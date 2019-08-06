import pygmsh
import meshio
import numpy as np


def mesh_el(eltype, points, lcar=0.1, filename='element.mesh'):
    """
    This function creates a GMSH geometry based on the input provided.
    In case 3 points are provided a triangle is created with a triangular 
    mesh. In case 4 points are provided a tetrahedron is created with a 
    tetrahedral mesh. The size of the mesh is set by the characteristic 
    length parameter lcar.
             
    Parameters
    ----------
    eltype:   string
              'CTRIA3' for triangular surface elements,
              'CQUAD4' for quadrilateral surface elements,
              'CTETRA' for tetrahedral volume elements,
              'CHEXA8' for hexahedronal volume elements.
    points:   array like
              Array containing the coordinates of the element
    lcar:     float
              Characteristic length value provided to gmsh for mesh sizing.
     
    filename: string
              file name of the mesh file in medit format which will be written    
         
    Returns 
    --------
    : array like
        Array with boolean entries stating True for those items
        which are considered to be part of the far field boundary
    
    """
    if eltype == "CTRIA3":
        sdim = 2
        pnt0, pnt1, pnt2 = points[0], points[1], points[2]
        geom = pygmsh.built_in.Geometry()
         
        p0 = geom.add_point(pnt0, lcar)
        p1 = geom.add_point(pnt1, lcar)
        p2 = geom.add_point(pnt2, lcar)
         
        l0 = geom.add_line(p0, p1)
        l1 = geom.add_line(p1, p2)
        l2 = geom.add_line(p2, p0)
    
        ll0 = geom.add_line_loop([l0, l1, l2])
        s0 = geom.add_plane_surface(ll0)

    if eltype == "CQUAD4":
        sdim = 2
        pnt0, pnt1, pnt2, pnt3 = points[0], points[1], points[2], points[3]
        geom = pygmsh.built_in.Geometry()

        p0 = geom.add_point(pnt0, lcar)
        p1 = geom.add_point(pnt1, lcar)
        p2 = geom.add_point(pnt2, lcar)
        p3 = geom.add_point(pnt3, lcar)
         
        l0 = geom.add_line(p0, p1)
        l1 = geom.add_line(p1, p2)
        l2 = geom.add_line(p2, p3)
        l3 = geom.add_line(p3, p0)
            
        ll0 = geom.add_line_loop([l0, l1, l2, l3])
        s0 = geom.add_plane_surface(ll0)    
    
    if eltype == "CTETRA":
        sdim = 3
        pnt0, pnt1, pnt2, pnt3 = points[0], points[1], points[2], points[3]

        geom = pygmsh.built_in.Geometry()

        p0 = geom.add_point(pnt0, lcar)
        p1 = geom.add_point(pnt1, lcar)
        p2 = geom.add_point(pnt2, lcar)
        p3 = geom.add_point(pnt3, lcar)
        
        l0 = geom.add_line(p0, p1)
        l1 = geom.add_line(p1, p2)
        l2 = geom.add_line(p2, p0)   
        l3 = geom.add_line(p0, p3)
        l4 = geom.add_line(p1, p3)
        l5 = geom.add_line(p2, p3)
        
        ll0 = geom.add_line_loop([l0, l1, l2])
        ll1 = geom.add_line_loop([l3, -l4, -l0])
        ll2 = geom.add_line_loop([l4, -l5, -l1])
        ll3 = geom.add_line_loop([l5, -l3, -l2])
       
        s0 = geom.add_plane_surface(ll0)
        s1 = geom.add_plane_surface(ll1)
        s2 = geom.add_plane_surface(ll2)
        s3 = geom.add_plane_surface(ll3)
        
        sl0 = geom.add_surface_loop([s0, s1, s2, s3])
        geom.add_volume(sl0)
        
    if eltype == "CHEXA8":
        sdim = 3
        pnt0, pnt1, pnt2, pnt3 = points[0], points[1], points[2], points[3]
        pnt4, pnt5, pnt6, pnt7 = points[4], points[5], points[6], points[7]
        geom = pygmsh.built_in.Geometry()
        
        p0 = geom.add_point(pnt0, lcar)
        p1 = geom.add_point(pnt1, lcar)
        p2 = geom.add_point(pnt2, lcar)
        p3 = geom.add_point(pnt3, lcar)
        p4 = geom.add_point(pnt4, lcar)
        p5 = geom.add_point(pnt5, lcar)
        p6 = geom.add_point(pnt6, lcar)
        p7 = geom.add_point(pnt7, lcar)
        
        l0 = geom.add_line(p0, p1)
        l1 = geom.add_line(p1, p2)
        l2 = geom.add_line(p2, p3)   
        l3 = geom.add_line(p3, p0)
        l4 = geom.add_line(p4, p5)
        l5 = geom.add_line(p5, p6)
        l6 = geom.add_line(p6, p7)
        l7 = geom.add_line(p7, p4)
        l8 = geom.add_line(p0, p4)   
        l9 = geom.add_line(p1, p5)
        l10 = geom.add_line(p2, p6)
        l11 = geom.add_line(p3, p7)
        
        ll0 = geom.add_line_loop([l0, l1, l2, l3])
        ll1 = geom.add_line_loop([-l7, -l6, -l5, -l4])
        ll2 = geom.add_line_loop([l8, l4, -l9, -l0])
        ll3 = geom.add_line_loop([l9, l5, -l10, -l1])
        ll4 = geom.add_line_loop([l10, l6, -l11, -l2])
        ll5 = geom.add_line_loop([-l3, l11, l7, -l8])
       
        s0 = geom.add_plane_surface(ll0)
        s1 = geom.add_plane_surface(ll1)
        s2 = geom.add_plane_surface(ll2)
        s3 = geom.add_plane_surface(ll3)
        s4 = geom.add_plane_surface(ll4)
        s5 = geom.add_plane_surface(ll5)
        
        sl0 = geom.add_surface_loop([s0, s1, s2, s3, s4, s5])
        geom.add_volume(sl0)
    
    'Generate mesh'
    points, cells, _, _, _ = pygmsh.generate_mesh(geom, 
                                                  optimize=False, num_lloyd_steps=0, 
                                                  dim=sdim, 
                                                  verbose=True)
 


    'Remove 3rd dimension in case of 2D element:'
    if (eltype == "CTRIA3") or (eltype == "CQUAD4"):
        points = np.delete(points, (2), axis=1)
        print "Meshed %i elements" % len(cells['triangle'] + 1)
    'Remove all keys in cells dict except for tetra information if sdim=3:'
    if (eltype == "CTETRA") or (eltype == "CHEXA8"):
        cells = {'tetra': cells['tetra']}
        print "Meshed %i elements" % (len(cells['tetra']) + 1)
    'Write resulting mesh to file'
    meshio.write(filename, points, cells, file_format='medit')
