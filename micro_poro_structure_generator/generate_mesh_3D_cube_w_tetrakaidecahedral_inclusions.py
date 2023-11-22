################################################################################
###                                                                          ###
### Created by Mahdi Manoochehrtayebi, 2020-2023                             ###
###                                                                          ###
### École Polytechnique, Palaiseau, France                                   ###
###                                                                          ###
################################################################################

import gmsh

import micro_poro_structure_generator as gen

################################################################################

def generate_mesh_3D_cube_w_tetrakaidecahedral_inclusions(
        mesh_filename,
        l,
        h,
        lcar):
    """Warning! This does not work…"""
    
    gmsh.initialize()
    gmsh.clear()
    model = gmsh.model
    occ = model.occ
    occ.mesh.ToleranceInitialDelaunay=1e-12
    model.mesh.ToleranceInitialDelaunay=1e-12

    # occ.mesh.MinimumCircleNodes=16
    # occ.mesh.MinimumCurveNodes=16
    # occ.mesh.MinimumCirclePoints=16
    # occ.mesh.MinimumCurvePoints=16
    # occ.mesh.MeshSizeFromCurvature=16
    # occ.mesh.MeshSizeFromPoints=0
    # occ.mesh.MeshSizeFromParametricPoints=0

    p1 = occ.addPoint(0, 0, l, lcar)
    p2 = occ.addPoint(0, 0, -l, lcar)
    p3 = occ.addPoint(l, 0, 0, lcar)
    p4 = occ.addPoint(0, -l, 0, lcar)
    p5 = occ.addPoint(-l, 0, 0, lcar)
    p6 = occ.addPoint(0, l, 0, lcar)

    s1 = occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(p1, p3), occ.addLine(p3, p6), occ.addLine(p6, p1)])])
    s2 = occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(p1, p6), occ.addLine(p6, p5), occ.addLine(p5, p1)])])
    s3 = occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(p1, p5), occ.addLine(p5, p4), occ.addLine(p4, p1)])])
    s4 = occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(p1, p4), occ.addLine(p4, p3), occ.addLine(p3, p1)])])
    s5 = occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(p3, p6), occ.addLine(p3, p2), occ.addLine(p2, p6)])])
    s6 = occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(p5, p6), occ.addLine(p6, p3), occ.addLine(p3, p5)])])
    s7 = occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(p4, p5), occ.addLine(p5, p2), occ.addLine(p2, p4)])])
    s8 = occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(p4, p3), occ.addLine(p3, p2), occ.addLine(p2, p4)])])

    occ.addVolume([occ.addSurfaceLoop([s1, s2, s3, s4, s5, s6, s7, s8])], 100)
    
    # setPeriodic(0)
    # setPeriodic(1)
    occ.mesh.ToleranceInitialDelaunay=1e-12
    occ.synchronize()
    gmsh.model.mesh.ToleranceInitialDelaunay=1e-12
    gmsh.fltk.run()
    model.mesh.generate()

    gmsh.write(mesh_filename + '.msh')

    gen.convert_msh_to_xml(mesh_filename)

    gmsh.write(mesh_filename+".vtk")
    gen.convert_vtk_to_xdmf(mesh_filename, dim=3)

    gmsh.finalize()
