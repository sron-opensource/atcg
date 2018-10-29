'''
Created on 25 dec. 2017
@author: henkvw

This Python 2.7 script is intended to create a thermal lumped element 
representation of an equivalent physical part represented by a FEM
mesh. This mesh can build from 2D triangles, 2D Quadrilateral elements
, 3D tetrahedral elements and 3D brick elements. In its current form 
the mesh needs to be provided in MSC Nastran format using the .nas ascii 
format.

The far-field method as implemented in this script is based on the following 
articles:

Automatic Conductor Generation for Thermal Lumped Parameter Models
by S. Appel, R. Patricio, H.P. de Koning and O. Pin
34th International Conference on Environmental Systems (ICES)
July 19-22, 2004, Colorado Springs, Colorado

and

Lumped Parameter Thermal Conductor Generation for 3D Geometry
by J.H. Strutt, N.J. Stock and C.J. Kirtley 
44th International Conference on Environmental Systems ICES-2014-297
13-17 July 2014, Tucson, Arizona

The script relies on the following existing software and Python
libraries:
gmsh:    Standalone mesh program used to create the meshes for the
         individual elements of the provided geometry.
pygmsh:  Python routines for automatic meshing of geometries in gmsh
         using python commands.
meshio:  Python module which is used to write the pygmsh mesh in a format
         which is recognized by sfepy.
sfepy:   The Python module which is used for implementing the required FEM
         routines. 
numpy:   Python module used for numerical routines and functions.
'''

import numpy as np
import os
import ntpath

from pyNastran.bdf.bdf import BDF

from __builtin__ import str
from lib.farfield import FarField
from lib.mesh import mesh_el

import datetime



def _rotate_2D(points):
    """
    This helper function transforms a triangle or quadrilateral element
    in 3D given by the coordinates p1, p2, p3,.. to an equivalent 
    triangle or quadrilateral in the x-y plane with vertices 
    p1n, p2n, p3n,.. These are returned in the form of an array.
   
    The method used is described on:
    https://en.wikipedia.org/wiki/Rotation_matrix

    Todo: Instead of taking the first 3 points and transforming the plane passing through them, a 
    plane should be fit to the points and then the plane (along with points) should be transformed 
    to become the x-y plane.
    
    Parameters
    ----------
    points: array like
            Array containing the entries of the points to be transformed
            using the form [[p1x, p1y, p1z],[p2x, p2y, p2z],[p3x, p3y, p3z]]
            for a triangle and for a quadrilateral element:
            [[p1x, p1y, p1z],[p2x, p2y, p2z],[p3x, p3y, p3z],[p4x, p4y, p4z]]
         
    Returns
    -------
    : array like
            Array with the transformed coordinates using the same format as the
            input array points
    """   
    if (len(points) == 3):
        pnt0, pnt1, pnt2 = np.array(points[0]), np.array(points[1]), np.array(points[2])
        '''
        Check if triangle is already aligned with xy plane and if it is return
        points with z set to 0:
        '''
        if (pnt2[2] == pnt0[2]) and (pnt1[2] == pnt0[2]):
            p0n = [pnt0[0], pnt0[1], 0.]
            p1n = [pnt1[0], pnt1[1], 0.]
            p2n = [pnt2[0], pnt2[1], 0.]
            return [p0n, p1n, p2n]
    else:
        pnt0, pnt1, pnt2, pnt3 = np.array(points[0]), np.array(points[1]), np.array(points[2]), np.array(points[3])
        if (pnt3[2] == pnt0[2]) and (pnt2[2] == pnt0[2]) and (pnt1[2] == pnt0[2]):
            p0n = [pnt0[0], pnt0[1], 0.]
            p1n = [pnt1[0], pnt1[1], 0.]
            p2n = [pnt2[0], pnt2[1], 0.]
            p3n = [pnt3[0], pnt3[1], 0.]
            return [p0n, p1n, p2n, p3n]
  
    '''
    If the element is not aligned create transformation matrices and
    perform transformation.
    '''
    v1 = pnt1 - pnt0
    v2 = pnt2 - pnt1
    n1 = np.cross(v1, v2)
    n2 = np.array([0., 0., 1.])
    U = np.cross(n1, n2) / np.linalg.norm(np.cross(n1, n2))
    C1 = np.dot(n1, n2) / (np.dot(np.linalg.norm(n1), np.linalg.norm(n2)))
    C2 = np.linalg.norm(np.cross(n1, n2)) / (np.dot(np.linalg.norm(n1), np.linalg.norm(n2)))
    M = np.array([[0, -U[2], U[1]],
                  [U[2], 0, -U[0]],
                  [-U[1], U[0], 0]])
    'Create rotation matrix from coefficients above'
    Rot = C1 * np.eye(3) + C2 * M + (1 - C1) * np.outer(U, np.transpose(U))
    'Apply found rotation to all points'
    points2D = []
    for pnt in points:
        'Replace z value by zero'
        pnt2D = Rot.dot(pnt)[0:2]
        pnt2D = np.append(pnt2D, 0.)
        'Add 2D point to result'
        points2D.append(pnt2D)
    
    return points2D


def _all_on_edge_3d(sdim, points, vertices):
    if (sdim == 2):
        return False
    p1, p2, p3, p4 = points[0], points[1], points[2], points[3]
    edges = [p2 - p1, p3 - p2, p1 - p3, p4 - p1, p4 - p2, p4 - p3]
    origin = [p1, p2, p3, p1, p2, p3]
    result = [False] * len(vertices)
    for i in range(len(vertices)):
        vertex = vertices[i]
        temp = [np.linalg.norm(np.cross(edges[j], (vertex - origin[j]))) for j in range(len(edges))]
        result[i] = not(np.all(temp))
    if np.all(result):
        return True
    else:
        return False



def _getIF(self):
    self.IFdict = {}
    for spc in self.model.spcs:
        for node in self.model.SPC(spc):
            s1 = str(node).split("SPC")[1]
            nid = int(s1.split()[1])
            ngid = s1.split()[3].rstrip('.')
            if ngid in self.IFdict:
                self.IFdict[ngid]['nodes'].append(nid)
            else:
                self.IFdict[ngid] = {'nodes': [nid], 'elements': {}}
    for key in self.IFdict:
        for nid in self.IFdict[key]['nodes']:
            tmp = self.model.get_node_id_to_element_ids_map()[nid]
            for i in tmp:
                inodes = [node.nid for node in self.model.Element(i).nodes]
                if np.sum([inodes.count(x) for x in self.IFdict[key]['nodes']]) > 1:
                    if i not in self.IFdict[key]['elements']:
                        self.IFdict[key]['elements'].setdefault(i, [])
        for eid in self.IFdict[key]['elements']:
            el = self.model.Element(eid)
            elnodes = [node.nid for node in el.nodes]
            if el.type == 'CTRIA3':
                sides = [[0, 1], [1, 2], [2, 0]]
            if el.type == 'CQUAD4':
                sides = [[0, 1], [1, 2], [2, 3], [3, 0]]
            for i in range(len(sides)):
                v1 = elnodes[sides[i][0]]
                v2 = elnodes[sides[i][1]]
                if (v1 in self.IFdict[key]['nodes']) and (v2 in self.IFdict[key]['nodes']):
                    self.IFdict[key]['elements'][eid].append([[v1, v2], 0, 0])      
    IFlist = sorted([int(i) for i in list(self.IFdict.keys())])  
    print '-------------------------------'
    print "Found %i interface groups:" % len(self.IFdict)
    print '-------------------------------'
    for i in IFlist:
        print "%s: %s" % (str(i), self.IFdict[str(i)])   


def _interfacesection(self, IFnamedict=None):
    """
    """
    section = "   /* Interfaces section: */\n"
    space = 80
    nIF = len(self.IFdict)
    y_off = (nIF - 1) / 2.
    n = 0
    IFlist = sorted([int(i) for i in list(self.IFdict.keys())])
    named = False
    if IFnamedict is not None:
        named = True
    for i in IFlist:
        if named:
            IFkey = str(i)
            IFname = IFnamedict[IFkey]['name']
            section += "   ThermalCryogenics.Interfaces.CryoHeatPortNamed %s " % IFname
            section += "annotation(\n"
            section += "      Placement(visible = true,\n" 
            section += "         transformation(\n"
            section += "            origin = {%3.1f, %3.1f},\n" % (IFnamedict[IFkey]['transformation']['origin'][0], IFnamedict[IFkey]['transformation']['origin'][1]) 
            section += "            extent = {{-10, -10}, {10, 10}},\n" 
            section += "            rotation = 0),\n"
            section += "      iconTransformation(\n"
            section += "         origin = {%3.1f, %3.1f},\n" % (IFnamedict[IFkey]['iconTransformation']['origin'][0], IFnamedict[IFkey]['iconTransformation']['origin'][1]) 
            section += "         extent = {{-10, -10}, {10, 10}},\n" 
            section += "         rotation = 0)));\n"
        else:
            section += "   ThermalCryogenics.Interfaces.CryoHeatPortNamed IF_%s " % str(i)
            section += "annotation(\n"
            section += "      Placement(visible = true,\n" 
            section += "         transformation(\n"
            section += "            origin = {0, %3.3f},\n" % (space * (y_off - n)) 
            section += "            extent = {{-10, -10}, {10, 10}},\n" 
            section += "            rotation = 0),\n"
            section += "      iconTransformation(\n"
            section += "         origin = {0, %3.3f},\n" % (space * (y_off - n)) 
            section += "         extent = {{-10, -10}, {10, 10}},\n" 
            section += "         rotation = 0)));\n"
        n += 1 
    return section


def _modelsection(self):
    """
    """
    section = "\nprotected\n"
    section += "\n   /* Heat capacitors: */\n"
    for (eid, el) in self.conv_el:
        section += "   ThermalCryogenics.Components.HeatCapacitors.HeatCapacitor %s_C_%i (" % (self.modelname, eid)
        section += "T(fixed = true, start = Tstart),\n"
        section += "   redeclare replaceable package Material = Material, "
        section += "v (displayUnit = \"m3\") = %e);\n" % self.vol[eid - 1]
    section += "\n   /* Heat conductances internal to part: */\n"
    for (eid1, el) in self.conv_el:
        for eid2 in range(eid1, self.maxelid + 1):
            if self.con_mat[eid1 - 1][eid2 - 1] != 0:
                print "Calculating F between %i and %i" % (eid1, eid2)
                F12 = self.F_mat[eid1 - 1][eid2 - 1]
                print "F12 : %e" % F12
                F21 = self.F_mat[eid2 - 1][eid1 - 1]
                print "F21 : %e" % F21
                F = 1 / (1 / F12 + 1 / F21)
                if self.PlatingT is not None:
                    print "Calculating Fpl between %i and %i" % (eid1, eid2)
                    Fpl12 = self.Fpl_mat[eid1 - 1][eid2 - 1]
                    print "Fpl12 : %e" % Fpl12
                    Fpl21 = self.Fpl_mat[eid2 - 1][eid1 - 1]
                    print "Fpl21 : %e" % Fpl21
                    Fpl = 1 / (1 / Fpl12 + 1 / Fpl21)
                section += "   ThermalCryogenics.Components.BulkThermalResistors.BulkThermalResistor %s_G_%i_%i (\n" % (self.modelname, eid1, eid2)
                section += "   redeclare replaceable package Material = Material, "
                section += "A (displayUnit = \"m2\") = %.2e, l = 1.0);\n" % F
                if self.PlatingT is not None:
                    section += "   ThermalCryogenics.Components.BulkThermalResistors.BulkThermalResistor %s_Gpl_%i_%i (\n" % (self.modelname, eid1, eid2)
                    section += "   redeclare replaceable package Material = Plating, "
                    section += "A (displayUnit = \"m2\") = %.2e, l = 1.0);\n" % Fpl
    section += "\n   /* Heat conductances from part to interfaces: */\n"
    for id, IF in self.IFdict.items():
        for eid in IF['elements']:
            for edge in IF['elements'][eid]:
                section += "   ThermalCryogenics.Components.BulkThermalResistors.BulkThermalResistor %s_G_IF%i_%i_%i (\n" % (self.modelname, eid, edge[0][0], edge[0][1])
                section += "   redeclare replaceable package Material = Material, "
                section += "A (displayUnit = \"m2\") = %.2e, l = 1.0);\n" % edge[1]
                if self.PlatingT is not None:
                    section += "   ThermalCryogenics.Components.BulkThermalResistors.BulkThermalResistor %s_Gpl_IF%i_%i_%i (\n" % (self.modelname, eid, edge[0][0], edge[0][1])
                    section += "   redeclare replaceable package Material = Plating, "
                    section += "A (displayUnit = \"m2\") = %.2e, l = 1.0);\n" % edge[2]
    return section


def _connectsection(self, IFnamedict=None):
    """
    """
    section = "   /* Connections section: */\n"
    section += "equation\n"
    section += "   /* Internal part connections: */\n"
    for (eid1, el) in self.conv_el:
        for eid2 in range(eid1, self.maxelid + 1):
            if self.con_mat[eid1 - 1][eid2 - 1] != 0:
                section += "   connect(%s_C_%i.port, %s_G_%i_%i.port_a);\n" % (self.modelname, eid1, self.modelname, eid1, eid2)
                section += "   connect(%s_G_%i_%i.port_b, %s_C_%i.port);\n" % (self.modelname, eid1, eid2, self.modelname, eid2)
                if self.PlatingT is not None:
                    section += "   connect(%s_C_%i.port, %s_Gpl_%i_%i.port_a);\n" % (self.modelname, eid1, self.modelname, eid1, eid2)
                    section += "   connect(%s_Gpl_%i_%i.port_b, %s_C_%i.port);\n" % (self.modelname, eid1, eid2, self.modelname, eid2)   
    section += "   /* Part connections to external interfaces: */\n"
    for id, IF in self.IFdict.items():
        for eid in IF['elements']:
            for edge in IF['elements'][eid]:
                section += "   connect(%s_C_%i.port, %s_G_IF%i_%i_%i.port_a);\n" % (self.modelname, eid, self.modelname, eid, edge[0][0], edge[0][1])
                if self.PlatingT is not None:
                    section += "   connect(%s_C_%i.port, %s_Gpl_IF%i_%i_%i.port_a);\n" % (self.modelname, eid, self.modelname, eid, edge[0][0], edge[0][1])   
                if IFnamedict is None:
                    section += "   connect(%s_G_IF%i_%i_%i.port_b, IF_%s);\n" % (self.modelname, eid, edge[0][0], edge[0][1], id)
                    if self.PlatingT is not None:
                        section += "   connect(%s_Gpl_IF%i_%i_%i.port_b, IF_%s);\n" % (self.modelname, eid, edge[0][0], edge[0][1], id)
                else:
                    IFkey = id
                    IFname = IFnamedict[IFkey]['name']
                    section += "   connect(%s_G_IF%i_%i_%i.port_b, %s);\n" % (self.modelname, eid, edge[0][0], edge[0][1], IFname)
                    if self.PlatingT is not None:
                        section += "   connect(%s_Gpl_IF%i_%i_%i.port_b, %s);\n" % (self.modelname, eid, edge[0][0], edge[0][1], IFname)
    return section
    

class FEMModel:
    """
    Class representing a FEM model which is based on a Nastran bdf file
    and which needs to be converted to a Modelica lumped element model.
    From the Nastran bdf file only CTRIA3 triangular plate elements, CQUAD4
    quadrilateral plate elements and CTETRA four-Sided solid Element are
    imported. The corresponding Modelica lumped element models are CTRIA3.mo,
    CQUAD4.mo and CTETRA.mo  
    """

    count = 0   # How many FEMModel objects are there?
    
    def __init__(self, NastranFile, meshfile):
        """
        Extend `FEMModel.__init__()` to append required input and output data.
        Keyword arguments:
        NastranFile -- Nastran input filename
        """

        self.NastranFile = NastranFile
        self.model = None
        self.conv_el = []
        self.maxelid = None
        self.NCTRIA3 = None
        self.NCQUAD4 = None
        self.NCTETRA = None
        self.NCHEXA8 = None
        self.IFdict = {}

        self._readNas()
        self._getIF()
        
        self.con_mat = np.zeros((self.maxelid, self.maxelid))
        self.vol = np.zeros(self.maxelid)
        self.F_mat = np.zeros((self.maxelid, self.maxelid))
        self.meshfile = meshfile
        self.calc_vol()
        self.PlatingT = None
        print "\nTotal volume of nastran file: %3.3e m^3" % np.sum(self.vol)

        FEMModel.count += 1

    def _readNas(self):
        """ Read the nastran file and extract the elements.
        """
            
        'Create the BDF object'
        print "Start reading Nastran file:"
        print "%s" % os.path.abspath(self.NastranFile) 
        self.model = BDF(debug=False)
        self.model.read_bdf(self.NastranFile, xref=True)
        self.conv_el = [(eid, el) for (eid, el) in sorted(self.model.elements.items()) 
                        if el.type in ['CTRIA3', 'CQUAD4', 'CTETRA', 'CHEXA8']]
        self.maxelid = np.max(self.model.element_ids)
        
        self.NCTRIA3 = len([(eid, el) for (eid, el) in self.conv_el if el.type == 'CTRIA3'])
        self.NCQUAD4 = len([(eid, el) for (eid, el) in self.conv_el if el.type == 'CQUAD4'])
        self.NCTETRA = len([(eid, el) for (eid, el) in self.conv_el if el.type == 'CTETRA'])
        self.NCHEXA8 = len([(eid, el) for (eid, el) in self.conv_el if el.type == 'CHEXA8'])
        
        print '-------------------------------'
        print '%-4i elements found in file which can be converted:' % len(self.conv_el)
        print '-------------------------------'
        print "%-4i CTRIA3 elements" % self.NCTRIA3
        print "%-4i CQUAD4 elements" % self.NCQUAD4
        print "%-4i CTETRA elements" % self.NCTETRA
        print "%-4i CHEXA8 elements" % self.NCHEXA8
        print '-------------------------------'
        print "Finished reading Nastran file."
        print '-------------------------------\n'


    def _getIF(self):
        """ Get interfaces.
        """
        for spc in self.model.spcs:
            for node in self.model.SPC(spc):
                s1 = str(node).split("SPC")[1]
                nid = int(s1.split()[1])
                ngid = s1.split()[3].rstrip('.')
                if ngid in self.IFdict:
                    self.IFdict[ngid]['nodes'].append(nid)
                else:
                    self.IFdict[ngid] = {'nodes': [nid], 'elements': {}}
        for key in self.IFdict:
            for nid in self.IFdict[key]['nodes']:
                tmp = self.model.get_node_id_to_element_ids_map()[nid]
                for i in tmp:
                    inodes = [node.nid for node in self.model.Element(i).nodes]
                    if np.sum([inodes.count(x) for x in self.IFdict[key]['nodes']]) > 1:
                        if i not in self.IFdict[key]['elements']:
                            self.IFdict[key]['elements'].setdefault(i, [])
            for eid in self.IFdict[key]['elements']:
                el = self.model.Element(eid)
                elnodes = [node.nid for node in el.nodes]
                if el.type == 'CTRIA3':
                    sides = [[0, 1], [1, 2], [2, 0]]
                if el.type == 'CQUAD4':
                    sides = [[0, 1], [1, 2], [2, 3], [3, 0]]
                for i in range(len(sides)):
                    v1 = elnodes[sides[i][0]]
                    v2 = elnodes[sides[i][1]]
                    if (v1 in self.IFdict[key]['nodes']) and (v2 in self.IFdict[key]['nodes']):
                        self.IFdict[key]['elements'][eid].append([[v1, v2], 0])      

        print self.IFdict

        
    def _connected(self, el_type_1, el_nodes_1, el_type_2, el_nodes_2):
        """
        """
        comb = {'CTRIA3CTRIA3': 2,
                'CTRIA3CQUAD4': 2,
                'CTRIA3CTETRA': 2,
                'CTRIA3CHEXA8': 2,
                'CQUAD4CTRIA3': 2,
                'CQUAD4CQUAD4': 2,
                'CQUAD4CTETRA': 2,
                'CQUAD4CHEXA8': 2,
                'CTETRACTRIA3': 2,
                'CTETRACQUAD4': 2,
                'CTETRACTETRA': 3,
                'CTETRACHEXA8': -1,
                'CHEXA8CTRIA3': 2,
                'CHEXA8CQUAD4': 2,
                'CHEXA8CTETRA': -1,
                'CHEXA8CHEXA8': 4,
                }
        mask = el_type_1 + el_type_2
        inters = len(np.intersect1d(el_nodes_1, el_nodes_2))
        if inters == comb[mask]:
            value = True
        else:
            value = False
        return value

    def _get_con(self):
        """
        """
        nid_to_eids_map = self.model.get_node_id_to_element_ids_map()
        for (eid, el) in sorted(self.model.elements.items()):
            con_el = []
            el_type_1 = el.type
            el_nodes_1 = el.node_ids
            con_el_node = []
            for elid in el.node_ids:
                [con_el_node.append(i) for i in nid_to_eids_map[elid]]      
            for eid2 in set(con_el_node):
                el_type_2 = self.model.Element(eid2).type
                el_nodes_2 = self.model.Element(eid2).node_ids
                if self._connected(el_type_1, el_nodes_1, el_type_2, el_nodes_2):
                    con_el.append(eid2)
            for con in con_el:
                self.con_mat[con - 1][eid - 1] = 1
                self.con_mat[eid - 1][con - 1] = 1
        self.ncon = np.count_nonzero(np.triu(self.con_mat))
        print "\nTotal number of internal connections: %i" % (np.count_nonzero(self.con_mat))
    
    def calc_vol(self):
        """
        """
        for (eid, el) in sorted(self.model.elements.items()):
            v, lcar = self._calc_el_vol(el, eid)
            self.vol[eid - 1] = v
             
    def _calc_el_vol(self, el, eid):
        """
        """
        if (el.type == 'CTRIA3') or (el.type == 'CQUAD4'):
            a = el.Area()
            t = el.Thickness()
            v = a * t * 1e-9
            lcar = np.sqrt((4. / 3.) * np.sqrt(3.) * (a * 1 / 115.))
        if el.type == 'CTETRA':
            v = el.Volume() * 1e-9
            lcar = (v * 3e6 * np.sqrt(2)) ** (1. / 3.)
        if el.type == 'CHEXA8':
            v = el.Volume() * 1e-9
            lcar = (v * 1.5e6 * np.sqrt(2)) ** (1. / 3.)
        self.vol[eid - 1] = v
        return (v, lcar)    

    def convert(self, verbose=False, PlatingT=None):
        """
        """
        if PlatingT is not None:
            self.PlatingT = PlatingT
            self.Fpl_mat = np.zeros((self.maxelid, self.maxelid))
            print "Adding data for plating thickness of %3.3e" % self.PlatingT
        self._get_con()
        nF = 0     
        for (eid1, el1) in self.conv_el:
            if eid1:
                el1 = self.model.elements[eid1]
                (vol, lcar) = self._calc_el_vol(el1, eid1)
                print "\nElement %i evaluated: type=%s, lcar=%f, vol=%e" % (eid1, el1.type, lcar, 
                                                                            vol)
                points = self.model.elements[eid1].get_node_positions()
                if (el1.type == "CTRIA3") or (el1.type == "CQUAD4"):
                    points = _rotate_2D(points)
                if (el1.type == "CTETRA") or (el1.type == "CHEXA8"):
                    el13d = True
                    t = None
                else:
                    el13d = False
                    t = el1.Thickness()
                print "Creating mesh..."
                mesh_el(el1.type, points, lcar=lcar, filename=self.meshfile)
                for eid2 in np.arange(1, self.maxelid + 1):
                    if eid2:
                        if self.con_mat[eid1 - 1][eid2 - 1] != 0:
                            n_ids_1 = self.model.elements[eid1].node_ids
                            n_ids_2 = self.model.elements[eid2].node_ids
                            el2 = self.model.elements[eid2]
                            boundary = np.in1d(n_ids_1, n_ids_2)
                            if (el2.type == "CTRIA3") or (el2.type == "CQUAD4"):
                                el22d = True
                            else:
                                el22d = False
                            if el13d and el22d:
                                print "3D to 2D connection detected, adjusting boundary."
                                nodes = []
                                for node in n_ids_1:
                                    if node not in n_ids_2:
                                        nodes.append(node)
                                Ftot = 0
                                for node in nodes:
                                    for i in range(len(boundary)):
                                        if n_ids_1[i] == node:
                                            boundary[i] = True
                                    try:
                                        F = FarField(el1.type, points, boundary, lcar, 0.005, 
                                                     self.meshfile, thickness=t, verbose=verbose)
                                        print "Handling element %i, node %i added, F=%f" % (eid2, node, F)
                                        Ftot += F
                                    except Exception as e: 
                                        print "Handling element %i, F determination failed!" % (eid2)
                                        print(e)
                                Favg = Ftot / len(nodes)
                                self.F_mat[eid1 - 1][eid2 - 1] = Favg
                                nF += 1
                                print "Handling element %i, done, Favg=%f" % (eid2, Favg)
                            try:
                                F = FarField(el1.type, points, boundary, lcar, 0.005, self.meshfile, thickness=t, verbose=verbose)
                                self.F_mat[eid1 - 1][eid2 - 1] = F
                                nF += 1
                                print "Step %i, handling element %i, F=%f" % (nF, eid2, F)
                                if self.PlatingT is not None:
                                    self.Fpl_mat[eid1 - 1][eid2 - 1] = (self.PlatingT / t) * F
                                    print "Step %i, handling plating element %i, Fpl=%f" % (nF, eid2, self.Fpl_mat[eid1 - 1][eid2 - 1])
                            except Exception as e: 
                                print "Handling element %i, F determination failed!" % (eid2)
                                print(e)
                for key in self.IFdict:
                    if eid1 in self.IFdict[key]['elements']:
                        IFbound = self.IFdict[key]['elements'][eid1]
                        for IF in IFbound:
                            edge = IF[0]
                            boundary = np.in1d(n_ids_1, edge)
                            try:
                                F = FarField(el1.type, points, boundary, lcar, 0.005, self.meshfile,
                                             thickness=t, verbose=verbose)
                                IF[1] = F
                                print "Handling element %i interface, F=%f" % (eid1, F)
                                if self.PlatingT is not None:
                                    IF[2] = (self.PlatingT / t) * F
                                    print "Plating element %i interface, Fpl=%f" % (eid1, IF[2])
                            except Exception as e: 
                                print "Handling element %i, interface failed!" % (eid1)
                                print(e)
                            
        self.vtotal = np.sum(self.vol)
        print '\nTotal converted volume: %e\n' % self.vtotal
    
    def ExportModelicaFile(self, FileName, WithinStr=None, MaterialStr=None, IFnamedict=None, DiagramStr=None, PlatingStr=None):
        """
        """
           
        print "Start writing Modelica file:"
        print "%s" % os.path.abspath(FileName)
        'Open the Modelica file buffer for writing'
        f = open(FileName, 'w')
        self.modelname = ntpath.basename(FileName)[:-3]
        
        'Write Modelica file header'
        if WithinStr is None: 
            header = "model %s\n" % self.modelname
        else:
            header = "within %s;\n\n" % WithinStr
            header += "model %s\n" % self.modelname
        header += "   /*\n"
        header += "   This Modelica file is created using the FarFieldConvert python script.\n"
        header += "   Nastran input file: %s\n" % os.path.abspath(self.NastranFile)
        header += "   Exported on %s, at %s\n" % (datetime.datetime.now().strftime("%d/%m/%Y"), datetime.datetime.now().strftime("%H:%M:%S"))
        header += "   The file Contains the following mesh elements converted into lumped elements:\n"
        header += "   %-5i CTRIA3 elements\n" % self.NCTRIA3
        header += "   %-5i CQUAD4 elements\n" % self.NCQUAD4
        header += "   %-5i CTETRA elements\n" % self.NCTETRA
        header += "   %-5i CHEXA8 elements\n\n" % self.NCHEXA8
        header += "   The following %i IF groups were converted from SPC temperature BC's:\n" % len(self.IFdict)
        
        sortedkeys = sorted(map(int, list(self.IFdict.keys())))
        for id in sortedkeys:
            key = str(id)
            if IFnamedict is None:
                header += "   %s: %s \n" % (key, str(self.IFdict[key]))
            else:
                header += "   %s: nodes: %s \n" % (IFnamedict[key]['name'], self.IFdict[key]['nodes'])
        if self.PlatingT is not None:
            header += "   A plating is added to the model with a thickness of %f mm\n" % self.PlatingT
            header += "   The heat capacity of this plating is neglected, while parallel conductances are added to the base material.\n"
        header += "   */\n"
        f.write(header)
        f.write("\n")
        
        'Create material declaration and start temperature statement:'
        Strmat = "   /* Material declaration: */\n"
        if MaterialStr is None:
            Strmat += "   replaceable package Material = ThermalCryogenics.Materials.Unity\n"
            Strmat += "   constrainedby Modelica.ThermalCryogenics.Interfaces.PartialMaterial;\n"
        else:
            if self.PlatingT is None:
                Strmat += "   replaceable package Material = %s\n" % MaterialStr
                Strmat += "   constrainedby Modelica.ThermalCryogenics.Interfaces.PartialMaterial;\n"
            else:
                Strmat += "   replaceable package Material = %s\n" % MaterialStr
                Strmat += "   constrainedby Modelica.ThermalCryogenics.Interfaces.PartialMaterial;\n"
                Strmat += "   replaceable package Plating = %s\n" % PlatingStr
                Strmat += "   constrainedby Modelica.ThermalCryogenics.Interfaces.PartialMaterial;\n"
        f.write(Strmat)
        StrTstart = "   parameter Modelica.SIunits.Temperature Tstart(displayUnit=\"K\") \"Start temperature of part\";\n"
        f.write(StrTstart)
        f.write("\n")
            
        'Create component statements in mo file'
        f.write('   /* Component section: */\n')
        IFStr = _interfacesection(self, IFnamedict)
        f.write(IFStr)
        ModelStr = _modelsection(self)
        f.write(ModelStr)
        
        'Create connection equations in mo file'
        f.write('\n')
        self.ConnectStr = _connectsection(self, IFnamedict)
        f.write(self.ConnectStr)
        
        'Add footer in mo file'
        if DiagramStr is None:
            l = len(self.IFdict) * 50.
            footer = "\n"
            footer += "annotation(Icon(coordinateSystem(initialScale = 0.1),\n"
            footer += "   graphics={\n"
            footer += "      Text(\n"
            footer += "         origin = {0, 0},\n" 
            footer += "         lineColor = {0, 0, 255},\n" 
            footer += "         extent = {{-150, 140}, {150, 100}},\n"
            footer += "         textString = \"%name\"),\n" 
            footer += "      Rectangle(\n"
            footer += "         origin = {0, 0},\n"
            footer += "         extent = {{-100, %.0f}, {100, %.0f}})}),\n" % (l, -l)
            footer += "   uses(Modelica(version = \"3.2.2\")));\n"

        else:
            footer = DiagramStr
        footer += "\nend %s;" % self.modelname 
        f.write(footer)
        'Close the mo file buffer'
        f.close()
        print "Finished"
