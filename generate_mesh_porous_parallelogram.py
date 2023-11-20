################################################################################
###                                                                          ###
### Created by Mahdi Manoochhertaybei, 2020-2023                             ###
###                                                                          ###
### École Polytechnique, Palaiseau, France                                   ###
###                                                                          ###
################################################################################

import math
import meshio
import gmsh
import os
import dolfin

################################################################################

def generate_mesh_porous_parallelogram(mesh_filename, width, angle, r, shift_x, shift_y):


    gmsh.initialize()
    gmsh.clear()
    model = gmsh.model
    occ = model.occ

        ############# Parameters #######################################################
    a = width
    R = r
    t = angle
    b = math.sin(t)*a
    c = a * math.cos(t)

    lcar = a/50
        ############# Functions ########################################################
    xmin = 0; xmax = c
    ymin = 0; ymax = b
    zmin = 0; zmax = 0
    e = 1e-6

    occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(occ.addPoint(0, 0, 0, lcar), occ.addPoint(a, 0, 0, lcar)),\
                                        occ.addLine(occ.addPoint(a, 0, 0, lcar), occ.addPoint(a+c, b, 0, lcar)),\
                                        occ.addLine(occ.addPoint(a+c, b, 0, lcar), occ.addPoint(c, b, 0, lcar)),\
                                        occ.addLine(occ.addPoint(c, b, 0, lcar), occ.addPoint(0, 0, 0, lcar))])], 1)

    # occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(occ.addPoint(0+3*r, 0+3*r, 0, lcar), occ.addPoint(a+3*r, 0+3*r, 0, lcar)),\
    #                                     occ.addLine(occ.addPoint(a+3*r, 0+3*r, 0, lcar), occ.addPoint(a+c+3*r, b+3*r, 0, lcar)),\
    #                                     occ.addLine(occ.addPoint(a+c+3*r, b+3*r, 0, lcar), occ.addPoint(c+3*r, b+3*r, 0, lcar)),\
    #                                     occ.addLine(occ.addPoint(c+3*r, b+3*r, 0, lcar), occ.addPoint(0+3*r, 0+3*r, 0, lcar))])], 1)
    
    occ.addDisk(xc=0+shift_x, yc=0+shift_y, zc=0, rx=R, ry=R, tag=2)
    occ.addDisk(xc=a+shift_x, yc=0+shift_y, zc=0, rx=R, ry=R, tag=3)
    occ.addDisk(xc=a+c+shift_x, yc=b+shift_y, zc=0, rx=R, ry=R, tag=4)
    occ.addDisk(xc=c+shift_x, yc=b+shift_y, zc=0, rx=R, ry=R, tag=5)
    occ.cut(objectDimTags=[(2, 1)], toolDimTags=[(2, 2), (2, 3), (2, 4), (2, 5)], tag=6)
    occ.synchronize()
    gmsh.model.addPhysicalGroup(2, [6])

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
                    gmsh.model.mesh.setPeriodic(2, [j[1]], [i[1]], [1, 0, 0, dx,\
                                                                    0, 1, 0, dy,\
                                                                    0, 0, 1, dz,\
                                                                    0, 0, 0, 1 ])

    gmsh.initialize()


    setPeriodic(0)
    setPeriodic(1)
    gmsh.model.mesh.setSize(dimTags=gmsh.model.getEntities(0), size=lcar)
    gmsh.model.mesh.generate(dim=2)



    # gmsh.model.mesh.generate()
    gmsh.write(mesh_filename + '.msh')
    os.system("gmsh -2 -o " + mesh_filename + ".msh -format msh22 " + mesh_filename + ".msh")
    os.system("dolfin-convert " + mesh_filename + ".msh " +  mesh_filename + ".xml")

    gmsh.write(mesh_filename+"-mesh.vtk")
    gmsh.finalize()

    mesh = meshio.read(mesh_filename+"-mesh.vtk")
    mesh.points = mesh.points[:, :2]
    meshio.write(mesh_filename+"-mesh.xdmf", mesh)

