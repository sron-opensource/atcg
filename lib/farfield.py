from sfepy.base.base import output, Struct
import sfepy.discrete as sfedis
from sfepy.discrete.conditions import Conditions, EssentialBC
from sfepy.terms import Term
from sfepy.solvers.ls import ScipyDirect
from sfepy.solvers.nls import Newton
from sfepy.postprocess.viewer import Viewer

import numpy as np


def FarField(eltype, points, boundary, lcar, epsilon, meshfile, thickness=None, verbose=False):
    """
    This function determines a geometric factor F within a single element. 
    The element type can be triangular, quadrilateral, tetrahedral or 
    hexahedral. For these types eltype is set to  "CTRIA3", "CQUAD4"
    "CTETRA" and "CHEXA8" respectively. The vertices of the element are 
    provided in the input parameter points. The input parameter boundary is of 
    boolean type and has the same size as points. Those points which are
    part of the conductive interface CI are flagged True. Epsilon sets
    a tolerance to determine which mesh vertices are considered part of
    the FF boundary. The meshing of the element is stored in a location
    provided by meshfile.
    
    Parameters
    ----------
    
    eltype:   string
              'CTRIA3' for triangular surface elements,
              'CQUAD4' for quadrilateral surface elements,
              'CTETRA' for tetrahedral volume elements.
              'CHEXA8' for hexahedral volume elements.
    points:   array like
              Array containing the coordinates of the element
    boundary: array like
              Array containing the entries of the boundary points
    lcar:     float
              Characteristic length value provided to gmsh for mesh sizing.
    epsilon:  float
              Numerical tolerance on criterion for far field
              boundary
    meshfile: string
              Filename and location for storing temporary mesh file
    verbose:  boolean
              Indicate whether intermediate results should be displayed
              or not
         
    Returns
    -------
    : array like
        Array with boolean entries stating True for those items
        on the boundary and False otherwise
    """  
    
    if verbose is True:
        output.set_output(quiet=False)
    else:
        output.set_output(quiet=True)
    
    if (eltype == "CTRIA3") or (eltype == "CQUAD4"):
        sdim = 2
    else:
        sdim = 3
        
    boundpnts = []
    for i in range(len(points)):
        if boundary[i]:
            boundpnts.append(points[i])
            
    mesh = sfedis.fem.Mesh.from_file(meshfile)
    
    domain = sfedis.fem.FEDomain('domain', mesh)
    
    c = sfedis.Material('c', val=1.0)
    
    omega = domain.create_region('Omega', 'all')
    
    if verbose is True:
        coors = mesh.coors
        fixed_vert = _is_on_bound(coors, bound=boundpnts, sdim=sdim, epsilon=epsilon)
        print "fixed vertices:"
        print fixed_vert
 
    
    
    is_on_bound = sfedis.Functions([sfedis.Function('_is_on_bound', _is_on_bound, 
                                    extra_args={'bound': boundpnts, 'sdim': sdim, 
                                                'epsilon': lcar / 100.}), ])
    fixed = domain.create_region('fixed', 'vertices by _is_on_bound', 'facet', 
                                 functions=is_on_bound, add_to_regions=True)
    

    field_t = sfedis.fem.Field.from_args('temperature', np.float64, 'scalar', omega, approx_order=2)
    t = sfedis.FieldVariable('t', 'unknown', field_t, 1)
    s = sfedis.FieldVariable('s', 'test', field_t, 1, primary_var_name='t')
     
    integral = sfedis.Integral('i', order=4)
     
    term1 = Term.new('dw_laplace(s, t)', integral, omega, s=s, t=t)
    term2 = Term.new('dw_volume_integrate(c.val, s)', integral, omega, c=c, s=s)   # heat source term for 1st step of far field
    eq = sfedis.Equation('temperature', term1 - term2)
    eqs = sfedis.Equations([eq])
     
    t_fixed = EssentialBC('t_fixed', fixed, {'t.0': 0.0})
     
    ls = ScipyDirect({})
    nls = Newton({'i_max': 1, 'eps_a': 1e-10}, lin_solver=ls)
     
    pb = sfedis.Problem('temperature', equations=eqs, nls=nls, ls=ls)
    pb.time_update(ebcs=Conditions([t_fixed, ]))
     
    temperature = pb.solve()
    out = temperature.create_output_dict()
    
    if verbose is True:
        pb.save_state('result.vtk', out=out)
        view = Viewer('result.vtk')
        view(is_wireframe=True, rel_scaling=1, is_scalar_bar=True)
        print "Maximum temperature: %f" % np.max(out['t'].data)
    
    
    data = [i[0] for i in out['t'].data]
     
    FF = _get_far(eltype, points, data, mesh, sdim, epsilon)
    str1 = ''.join(str(v) + ', ' for v in FF)[:-2]
    try:
        far = domain.create_region('far', 'vertex %s' % str1, 'facet', add_to_regions=True)
    except Exception as e: 
        print "Far field region creation failed!"
        print(e)
        t.reset()
        s.reset()
        return
     
    area_source = pb.evaluate('d_surface.3.far(t)')
    
    fluxval = 1.0 / (area_source)
    c2 = sfedis.Material('c2', val=fluxval)             # So that total heat at the far field is 1W equally distributed over all elements
     
    term1A = Term.new('dw_laplace(c.val, s, t)', integral, omega, c=c, s=s, t=t)
    term2A = Term.new('dw_surface_integrate(c2.val, s)', integral, far, c2=c2, s=s)
    eq2 = sfedis.Equation('temperature2', term1A - term2A)
    eqs2 = sfedis.Equations([eq2])
       
    pb2 = sfedis.Problem('temperature2', equations=eqs2, nls=nls, ls=ls)
    pb2.time_update(ebcs=Conditions([t_fixed, ]))
     
    temperature2 = pb2.solve()
    out2 = temperature2.create_output_dict()
    volume = pb2.evaluate('d_volume.3.Omega(t)')
 
    t_int = pb2.evaluate('ev_volume_integrate.3.Omega(t)')

    avg_t = t_int / volume
    F = 1.0 / avg_t
    
    if verbose is True:
            print "Average temperature: %f" % avg_t
    
    if thickness:
        'Correction factor 1e-3 is due to geometry in mm instead of m'
        F = F * thickness * 1e-3
    
    if verbose is True:
        pb.save_state('result.vtk', out=out2)
        view = Viewer('result.vtk')
        view(is_wireframe=True, rel_scaling=1, is_scalar_bar=True)
    
    t.reset()
    s.reset()
    
    return F



def _is_on_bound(coors, domain=None, bound=None, sdim=2, epsilon=1e-6):  # @UnusedVariable
    """
    Helper function to determine which of the provided points in coors
    are part of the boundary provided in bound. All points within a 
    range smaller than epsilon are flagged True in the returned result. 
    This array has a length equal to the length of coors and it can be 
    used in a numpy where statement to filter those vertices which are
    on the boundary bound. The parameter domain is needed to match the
    function template for the sfepy vertex by function command, however
    the parameter is not actually used within the function itself. 
    
    Parameters
    ----------
    
    coors:   array like
             Array containing the coordinates to be evaluated
    domain:  Sfepy domain object
             Default set to None, required for Sfepy function template
    bound:   array like
             Array containing the entries of the boundary points
    sdim:    Integer
             Space dimension of the domain: 2 or 3
    epsilon: float
             Numerical tolerance
         
    Returns
    -------
    : array like
        Array with boolean entries stating True for those items
        on the boundary and False otherwise
    """  
    
    cri = epsilon * np.ones(len(coors))
    
    

    idim = len(bound)
    if idim == 2:
        if sdim == 2:
            p1, p2 = bound[0][0:2], bound[1][0:2]
        else:
            p1, p2 = bound[0], bound[1]
        v1 = np.array(p2) - np.array(p1)
        res = [np.linalg.norm(np.cross(v1, (pnt - p1))) for pnt in coors]
        flag = np.nonzero(res < cri)[0]
    else:
        p1, p2, p3 = bound[0], bound[1], bound[2]
        v1 = np.array(p2) - np.array(p1)
        v2 = np.array(p3) - np.array(p1)
        n1 = np.cross(v1, v2)
        tmp = [np.cross(v1, (pnt - p1)) for pnt in coors]
        res = [np.linalg.norm(np.cross(n1, n2)) for n2 in tmp]
        flag = np.nonzero(res < cri)[0]
        
    return flag


def _get_far(eltype, points, data, mesh, sdim, epsilon=0.005):
    """
    Get the vertices which belong to the far boundary gamma, using the
    criterion t_gamma > (1-epsilon)*t_max. Each facet (edge for 2D,
    face for 3D) connected to these vertices is considered as part of
    the boundary gamma.
       
    Parameters
    ----------
    
    eltype:   string
             'CTRIA3' for triangular surface elements,
             'CQUAD4' for quadrilateral surface elements,
             'CTETRA' for tetrahedral volume elements.
             'CHEXA8' for hexahedral volume elements.
    points:  array like
             Array containing the coordinates of the element vertices
    data:    array like
             Array containing the temperature results from step 1
    mesh:    array like
             Sfepy mesh object
    sdim:    Integer
             Space dimension of the domain: 2 or 3
    epsilon: float
             Numerical tolerance
         
    Returns
    -------
    : array like
        Array with boolean entries stating True for those items
        which are considered to be part of the far field boundary
    """
    'Place all vertices which meet criterion in v_res:'
    cri = np.ones(len(data)) * (1 - epsilon) * np.max(data)
    v_res = np.where(data > cri)[0]
    


    'Get all mesh connections and cells in both 2D and 3D:'
    if sdim == 2:
        conn, cells = mesh.get_conn('2_3', ret_cells=True)
        comb = [[0, 1], [1, 2], [2, 0]]
    else:
        conn, cells = mesh.get_conn('3_4', ret_cells=True)
        comb = [[0, 1, 2], [0, 3, 1], [1, 3, 2], [2, 3, 0]]
    
    '''
    Step 1:
    Create list cell_sel with cells which have vertices that are 
    part of the far field boundary:
    '''
    cell_sel = []
    for vertex in v_res:
        vert = np.ones(np.shape(conn)) * vertex
        match = np.where(conn == vert)[0]
        for i in match:
            cell = cells[i]
            if cell not in cell_sel:
                cell_sel.append(cell)
    
    '''
    Important: This section needs re-writing
    Case which is considered is incomplete
    '''       
    if (sdim == 3) and (len(v_res) < 10):

        for vertex in v_res:
            coor = mesh.coors[vertex]
            is_on_el_vertex = np.any(3 == (0 == (points - coor)).sum(1))
            if is_on_el_vertex:
                v_res = [vertex]
 
    'Check for CQUAD4 element type'    
    isCQUAD4 = eltype == "CQUAD4"   
    
    'Check for single vertex in v_res'    
    onevertex = len(v_res) == 1

    'Check for sdim==3 if all vertes in v_res are on edges'
    vertices = [mesh.coors[i] for i in v_res]
    # allonedge3d = _all_on_edge_3d(sdim, points, vertices)
    allonedge3d = False
    if (onevertex or allonedge3d or isCQUAD4):
        '''
        Add all adjacent elements of dimension sdim-1 to the far field
        result which are at the boundary of the element.
        Check for each facet of the cells in cell_sel if this
        facet contains a vertex from v_res and is on the boundary
        of the element
        ''' 
        far_vrt = []
        for cell in cell_sel:
            if sdim == 2:
                for facet in comb:
                    k, l = facet[0], facet[1] 
                    facet_vert = [conn[cell][k], conn[cell][l]]
                    if set(facet_vert).intersection(set(v_res)):
                        facet_coors = [mesh.coors[i] for i in facet_vert]
                        if _is_on_outer_bounds(facet_coors, eltype, points, sdim):
                            far_vrt.append(facet_vert)    
            if sdim == 3:
                for facet in comb:
                    k, l, m = facet[0], facet[1], facet[2] 
                    facet_vert = [conn[cell][k], conn[cell][l], conn[cell][m]]
                    if set(facet_vert).intersection(set(v_res)):
                        facet_coors = [mesh.coors[i] for i in facet_vert]
                        if _is_on_outer_bounds(facet_coors, eltype, points, sdim):
                            far_vrt.append(facet_vert)
    
        far_set = set([val for sublist in far_vrt for val in sublist])
        flag = list(far_set)
        
    else:
        '''
        Important: This section is incomplete: above only test
        is for single vertex or all vertices on edge and sdim=3
        This needs to be further completed...
        '''
        flag = v_res
    return flag



def _is_on_outer_bounds(coors, eltype, points, sdim):
    """
    Helper function to determine which of the provided points in coors
    are part of the element boundaries determined by element points of
    an element of type eltype.
    All points within a range smaller than epsilon are flagged True 
    in the returned result. This array has a length equal to the length
    of coors and it can be used in a numpy where statement to filter 
    those vertices which are on the boundary bound. 
       
    Parameters
    ----------
    
    coors:   array like
             Array containing the coordinates to be evaluated
    eltype:   string
             'CTRIA3' for triangular surface elements,
             'CQUAD4' for quadrilateral surface elements,
             'CTETRA' for tetrahedral volume elements.
             'CHEXA8' for hexahedral volume elements.
    points:  array like
             Array containing the coordinates of the element vertices
    epsilon: float
             Numerical tolerance
         
    Returns
    -------
    : bool
        stating True if coors are on element boundary
    """

    on_outer_bound = False
    if sdim == 2:
        if eltype == 'CTRIA3':
            comb = [[0, 1], [1, 2], [2, 0]]
        elif eltype == 'CQUAD4':
            comb = [[0, 1], [1, 2], [2, 3], [3, 0]]
        else:
            print "unexpected 2D element"
            return
        for facet in comb:
            k, l = facet[0], facet[1] 
            bound = [points[k][:], points[l][:]]
            if len(_is_on_bound(coors, domain=None, bound=bound, sdim=sdim)) == sdim:
                on_outer_bound = True    
        
    if sdim == 3:
        if eltype == 'CTETRA':
            comb = [[0, 1, 2], [0, 3, 1], [1, 3, 2], [2, 3, 0]]
        elif eltype == 'CHEXA8':
            comb = [[0, 1, 2, 3], [4, 5, 6, 7], [0, 1, 4, 7], [1, 2, 5, 4], [2, 3, 6, 5], [3, 0, 7, 6]]
        else:
            print "unexpected 3D element"
            return
        for facet in comb:
            k, l, m = facet[0], facet[1], facet[2] 
            bound = [points[k][:], points[l][:], points[m][:]]
            if len(_is_on_bound(coors, domain=None, bound=bound, sdim=sdim)) == sdim:
                on_outer_bound = True
    
    return on_outer_bound
