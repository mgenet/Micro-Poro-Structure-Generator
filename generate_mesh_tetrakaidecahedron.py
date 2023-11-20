################################################################################
###                                                                          ###
### Created by Mahdi Manoochhertaybei, 2020-2023                             ###
###                                                                          ###
### École Polytechnique, Palaiseau, France                                   ###
###                                                                          ###
################################################################################

import meshio
import gmsh
import os
import dolfin

################################################################################

def generate_mesh_tetrakaidecahedron(mesh_filename, l, h):
    
    ###################################################################################################################################################################
    l =1
    h = 0.1
    lcar = h/5

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

    
    def setPeriodic(coord):
        # From https://gitlab.onelab.info/gmsh/gmsh/-/issues/744
        smin = gmsh.model.getEntitiesInBoundingBox(xmin - e, ymin - e, zmin - e,
                                                    (xmin + e) if (coord == 0) else (xmax + e),
                                                    (ymin + e) if (coord == 1) else (ymax + e),
                                                    (zmin + e) if (coord == 2) else (zmax + e),
                                                    2)
        dx = (xmax - xmin) if (coord == 0) else 0
        dy = (ymax - ymin) if (coord == 1) else 0
        dz = (zmax - zmin) if (coord == 2) else 0

        for i in smin:
            bb = gmsh.model.getBoundingBox(i[0], i[1])
            bbe = [bb[0] - e + dx, bb[1] - e + dy, bb[2] - e + dz,
                    bb[3] + e + dx, bb[4] + e + dy, bb[5] + e + dz]
            smax = gmsh.model.getEntitiesInBoundingBox(bbe[0], bbe[1], bbe[2],
                                                        bbe[3], bbe[4], bbe[5])
            for j in smax:
                bb2 = list(gmsh.model.getBoundingBox(j[0], j[1]))
                bb2[0] -= dx; bb2[1] -= dy; bb2[2] -= dz;
                bb2[3] -= dx; bb2[4] -= dy; bb2[5] -= dz;
                if ((abs(bb2[0] - bb[0]) < e) and (abs(bb2[1] - bb[1]) < e) and
                    (abs(bb2[2] - bb[2]) < e) and (abs(bb2[3] - bb[3]) < e) and
                    (abs(bb2[4] - bb[4]) < e) and (abs(bb2[5] - bb[5]) < e)):
                    gmsh.model.mesh.setPeriodic(1, [j[1]], [i[1]], [1, 0, 0, dx,\
                                                                    0, 1, 0, dy,\
                                                                    0, 0, 1, dz,\
                                                                    0, 0, 0, 1 ])

    # setPeriodic(0)
    # setPeriodic(1)
    occ.mesh.ToleranceInitialDelaunay=1e-12
    occ.synchronize()
    gmsh.model.mesh.ToleranceInitialDelaunay=1e-12
    gmsh.fltk.run()
    model.mesh.generate()

    gmsh.write(mesh_filename + '.msh')
    gmsh.write(mesh_filename+"-mesh.vtk")
    gmsh.finalize()

    os.system("gmsh -2 -o " + mesh_filename + ".msh -format msh22 " + mesh_filename + ".msh")
    os.system("dolfin-convert " + mesh_filename + ".msh " +  mesh_filename + ".xml")

    mesh = meshio.read(mesh_filename+"-mesh.vtk")
    mesh.points = mesh.points[:, :2]
    meshio.write(mesh_filename+"-mesh.xdmf", mesh)
