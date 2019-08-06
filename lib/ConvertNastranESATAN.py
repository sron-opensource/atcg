'''
Created on 11 okt. 2018

@author: henkvw
'''

from pyNastran.bdf.bdf import BDF
import os
import numpy as np

def readnas(NastranFile, ExpFileName, ModelName):
    """
    """
        
    'Create the BDF object'
    print "Start reading Nastran file:"
    print "%s" % os.path.abspath(NastranFile) 
    model = BDF(debug=False)
    model.read_bdf(NastranFile, xref=True)
    conv_el = [(eid,el) for (eid,el) in sorted(model.elements.items()) if el.type in ['CTRIA3', 'CQUAD4', 'CTETRA', 'CHEXA8']]
    maxelid = np.max(model.element_ids)
    
    NCTRIA3 = len([(eid,el) for (eid,el) in conv_el if el.type == 'CTRIA3'])
    NCQUAD4 = len([(eid,el) for (eid,el) in conv_el if el.type == 'CQUAD4'])
    NCTETRA = len([(eid,el) for (eid,el) in conv_el if el.type == 'CTETRA'])
    NCHEXA8 = len([(eid,el) for (eid,el) in conv_el if el.type == 'CHEXA8'])
    
    print '-------------------------------'
    print '%-4i elements found in file which can be converted:' % len(conv_el)
    print '-------------------------------'
    print "%-4i CTRIA3 elements" % NCTRIA3
    print "%-4i CQUAD4 elements" % NCQUAD4
    print "%-4i CTETRA elements" % NCTETRA
    print "%-4i CHEXA8 elements" % NCHEXA8
    print '-------------------------------'
    print "Finished reading Nastran file."
    print '-------------------------------\n'
    
    print "Start writing Modelica file:"
    print "%s" % os.path.abspath(ExpFileName)
    'Open the Modelica file buffer for writing'
    f = open(ExpFileName, 'w')
    header = 'BEGIN_MODEL %s\n\n' %ModelName
    header += 'BULK Bulk_Alu6061T651;\n'
    header += 'Bulk_Alu6061T651 = [1.000, 1.000, 1.000];\n'
    header += '\n\n'
    f.write(header)
    for eid, el in conv_el:
        elname = '%s_%s_%04d' %(ModelName, el.type, eid)
        points = model.elements[eid].get_node_positions()
        thickness =  model.elements[eid].Thickness()
        section = 'SHELL %s;\n' %elname
        if el.type == 'CQUAD4':
            section += '%s = SHELL_QUADRILATERAL (\n' %elname
            section += 'point1 = [%3.6f, %3.6f, %3.6f],\n' %(points[0][0]/1e3, points[0][1]/1e3, points[0][2]/1e3)
            section += 'point2 = [%3.6f, %3.6f, %3.6f],\n' %(points[1][0]/1e3, points[1][1]/1e3, points[1][2]/1e3)
            section += 'point3 = [%3.6f, %3.6f, %3.6f],\n' %(points[2][0]/1e3, points[2][1]/1e3, points[2][2]/1e3)
            section += 'point4 = [%3.6f, %3.6f, %3.6f],\n' %(points[3][0]/1e3, points[3][1]/1e3, points[3][2]/1e3)
            section += 'nbase1 = %i,\n' %(10000 + eid)
            section += 'nbase2 = %i,\n' %(10000 + eid)
            section += 'bulk = Bulk_Alu6061T651,\n'
            section += 'thick = %3.6f);\n\n' %thickness
        if el.type == 'CTRIA3':
            section += '%s = SHELL_TRIANGLE (\n' %elname
            section += 'point1 = [%3.6f, %3.6f, %3.6f],\n' %(points[0][0]/1e3, points[0][1]/1e3, points[0][2]/1e3)
            section += 'point2 = [%3.6f, %3.6f, %3.6f],\n' %(points[1][0]/1e3, points[1][1]/1e3, points[1][2]/1e3)
            section += 'point3 = [%3.6f, %3.6f, %3.6f],\n' %(points[2][0]/1e3, points[2][1]/1e3, points[2][2]/1e3)
            section += 'nbase1 = %i,\n' %(10000 + eid)
            section += 'nbase2 = %i,\n' %(10000 + eid)
            section += 'bulk = Bulk_Alu6061T651,\n'
            section += 'thick = %3.6f);\n\n' %(thickness/1e3)
        f.write(section)
        print eid
        print el.type
        points = model.elements[eid].get_node_positions()
        print points
        print model.elements[eid].Thickness()
    
    footer = 'END_MODEL\n'
    f.write(footer)
    f.close()
    
    
filename = 'c:/Users/henkvw/ownCloud/Athena/Thermal modelling/ForCSL/T2_332_G_5030_th_2.nas'
expfile = 'c:/Users/henkvw/ownCloud/Athena/Thermal modelling/ForCSL/T2_332_G_5030_th_2.erg'
ModelName = 'T2_332_G_5030_th_2'

readnas(filename, expfile, ModelName)