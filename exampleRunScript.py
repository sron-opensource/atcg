from fem_model import FEMModel
import warnings
import sys
import datetime

if __name__ == "__main__":

    "Start scripting file here"
    "----------------------------------------------------------------------"

    if not sys.warnoptions:
        warnings.simplefilter("ignore")

    meshfile = './FEM/test.mesh'

    epsilon = 0.005

    # eltype = "CTETRA"
    # eltype = "CTRIA3"
    eltype = "CQUAD4"
    # eltype = "CHEXA8"


    if eltype == "CTRIA3":
        p0 = [0., 0., 0.]
        p1 = [1., 0., 0.]
        p2 = [0., 1., 0.]
        points = [p0, p1, p2]
        bound = [True, True, False]

    if eltype == "CTETRA":
        p0 = [0., 0., 0.]
        p1 = [1., 0., 0.]
        p2 = [0., 1., 0.]
        p3 = [0., 0., 1.]
        points = [p0, p1, p2, p3]
        bound = [True, True, True, False]

    if eltype == "CQUAD4":
        p0 = [0., 0., 0.]
        p1 = [1., 0., 0.]
        p2 = [1.5, 1., 0.]
        p3 = [0., 1., 0.]
        points = [p0, p1, p2, p3]
        bound = [True, True, False, False]
        
    if eltype == "CHEXA8":
        p0 = [0., 0., 0.]
        p1 = [0.3, 0., 0.]
        p2 = [0.3, 0.2, 0.]
        p3 = [0., 0.2, 0.]
        p4 = [0., 0., 0.1]
        p5 = [0.3, 0., 0.1]
        p6 = [0.3, 0.2, 0.1]
        p7 = [0., 0.2, 0.1]
        points = [p0, p1, p2, p3, p4, p5, p6, p7]
        bound = [True, False, False, True, True, False, False, True]

    # mesh_el(eltype, points)
    # lcar = 0.1
    # mesh_el(eltype, points, lcar = 0.1, filename='./FEM/test.mesh')
    # G = FarField(eltype, points, bound, lcar, epsilon, meshfile, verbose = True)


    tstart = datetime.datetime.now()
    print "Started at: %s" % tstart
    # nas_filename = 'CTRIA3_4.nas'
    # nas_filename = r'CTRIA3_16.nas'
    # nas_filename = 'CQUAD4_2.nas'
    # nas_filename = r'CQUAD4_8.nas'
    # nas_filename = r'CQUAD4_50.nas'
    # nas_filename = 'Rect_4_tri.nas'
    nas_filename = 'Rect_8_quad.nas'
    # nas_filename = "Rect_20_irtri.nas"
    # nas_filename = "Rect_12_irquad.nas"

    import_dir = 'NastranImport/'
    export_dir = 'Modelica'
    Model = FEMModel(import_dir + nas_filename, 'element.mesh')
    mo_filename = nas_filename[:-4] + ".mo"
    Model.convert(verbose=False)
    Model.ExportModelicaFile(export_dir + mo_filename)
    # 
     
    tend = datetime.datetime.now()
    print "\nEnded at: %s" % tend
    duration = tend - tstart
    print "Duration: %s" % duration
