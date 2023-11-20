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
import pickle
import dolfin
import numpy
import scipy
import matplotlib.pyplot as plt

################################################################################

def generate_mesh_hexagonal_inclusion(mesh_filename, phi, R):
    
    gmsh.initialize()
    gmsh.clear()
    model = gmsh.model
    occ = model.occ

    ##########################################################
    h = R*(1 - phi)
    L = R + h
    lcar = h/10

    def hexagon_generator(x, y, k):
        occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(occ.addPoint(x, y + L, 0, lcar), occ.addPoint(x + L * math.sqrt(3)/2, y + L/2, 0, lcar)),\
                                               occ.addLine(occ.addPoint(x + L * math.sqrt(3)/2, y + L/2, 0, lcar), occ.addPoint(x + L * math.sqrt(3)/2, y +(-L/2), 0, lcar)),\
                                               occ.addLine(occ.addPoint(x+ L * math.sqrt(3)/2, y + (-L/2), 0, lcar), occ.addPoint(x, y+(-L), 0, lcar)),\
                                               occ.addLine(occ.addPoint(x, y +(-L), 0, lcar), occ.addPoint(x+(-L) * math.sqrt(3)/2, y +(-L/2), 0, lcar)),\
                                               occ.addLine(occ.addPoint(x+ (-L) * math.sqrt(3)/2, y + (-L/2), 0, lcar), occ.addPoint(x  -L * math.sqrt(3)/2, y + L/2, 0, lcar)),\
                                               occ.addLine(occ.addPoint(x -L * math.sqrt(3)/2, y + L/2, 0, lcar), occ.addPoint(x, y + L, 0, lcar))])], k)

        occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(occ.addPoint(x, y + R, 0, lcar), occ.addPoint(x + R * math.sqrt(3)/2, y + R/2, 0, lcar)),\
                                               occ.addLine(occ.addPoint(x + R * math.sqrt(3)/2, y + R/2, 0, lcar), occ.addPoint(x + R * math.sqrt(3)/2, y +(-R/2), 0, lcar)),\
                                               occ.addLine(occ.addPoint(x+ R * math.sqrt(3)/2, y + (-R/2), 0, lcar), occ.addPoint(x, y+(-R), 0, lcar)),\
                                               occ.addLine(occ.addPoint(x, y +(-R), 0, lcar), occ.addPoint(x+(-R) * math.sqrt(3)/2, y +(-R/2), 0, lcar)),\
                                               occ.addLine(occ.addPoint(x+ (-R) * math.sqrt(3)/2, y + (-R/2), 0, lcar), occ.addPoint(x  -R * math.sqrt(3)/2, y + R/2, 0, lcar)),\
                                               occ.addLine(occ.addPoint(x -R * math.sqrt(3)/2, y + R/2, 0, lcar), occ.addPoint(x, y + R, 0, lcar))])], k+1)


        occ.cut(objectDimTags=[(2, k)], toolDimTags=[(2, k+1)], tag=k+2)

    # occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(occ.addPoint(0, L, 0, lcar), occ.addPoint(L * math.sqrt(3)/2, L/2, 0, lcar)),\
    #                                        occ.addLine(occ.addPoint(L * math.sqrt(3)/2, L/2, 0, lcar), occ.addPoint(L * math.sqrt(3)/2, -L/2, 0, lcar)),\
    #                                        occ.addLine(occ.addPoint(L * math.sqrt(3)/2, -L/2, 0, lcar), occ.addPoint(0, -L, 0, lcar)),\
    #                                        occ.addLine(occ.addPoint(0, -L, 0, lcar), occ.addPoint(-L * math.sqrt(3)/2, -L/2, 0, lcar)),\
    #                                        occ.addLine(occ.addPoint(-L * math.sqrt(3)/2, -L/2, 0, lcar), occ.addPoint(-L * math.sqrt(3)/2, L/2, 0, lcar)),\
    #                                        occ.addLine(occ.addPoint(-L * math.sqrt(3)/2, L/2, 0, lcar), occ.addPoint(0, L, 0, lcar))])], 1)

    # occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(occ.addPoint(0, R, 0, lcar), occ.addPoint(R * math.sqrt(3)/2, R/2, 0, lcar)),\
    #                                        occ.addLine(occ.addPoint(R * math.sqrt(3)/2, R/2, 0, lcar), occ.addPoint(R * math.sqrt(3)/2, -R/2, 0, lcar)),\
    #                                        occ.addLine(occ.addPoint(R * math.sqrt(3)/2, -R/2, 0, lcar), occ.addPoint(0, -R, 0, lcar)),\
    #                                        occ.addLine(occ.addPoint(0, -R, 0, lcar), occ.addPoint(-R * math.sqrt(3)/2, -R/2, 0, lcar)),\
    #                                        occ.addLine(occ.addPoint(-R * math.sqrt(3)/2, -R/2, 0, lcar), occ.addPoint(-R * math.sqrt(3)/2, R/2, 0, lcar)),\
    #                                        occ.addLine(occ.addPoint(-R * math.sqrt(3)/2, R/2, 0, lcar), occ.addPoint(0, R, 0, lcar))])], 2)

    # occ.cut(objectDimTags=[(2, 1)], toolDimTags=[(2, 2)], tag=3)
    hexagon_generator(0, 0, 10)
    hexagon_generator((2*R + h)*math.sqrt(3)/2, 0, 20)
    hexagon_generator((2*R + h)*math.sqrt(3)/4, (2*R + h)*3/4, 30)
    hexagon_generator((2*R + h)*math.sqrt(3)/4, -(2*R + h)*3/4, 40)
    occ.fuse([(2,12)], [(2,22), (2, 32), (2, 42)], 50)
    shift = 0
    # shift = -R/5
    xmin = 0; dx = (2*R+h)*math.sqrt(3)/2; xmax = xmin + dx
    ymin = -L/2+shift; dy = (2*R+h)*3/2; ymax = ymin + dy
    zmin = 0; zmax = 0
    e = 1e-6
    occ.addRectangle(x=xmin, y=ymin-dy, z=0, dx=-dx, dy=2*dy, tag=100)
    occ.addRectangle(x=xmax, y=ymin-dy, z=0, dx=dx, dy=2*dy, tag=200)
    occ.addRectangle(x=xmin, y=ymin- R/2, z=0, dx=dx, dy=-dy, tag=300)
    occ.addRectangle(x=xmin, y=ymax-R/2, z=0, dx=dx, dy=dy, tag=400)
    occ.cut(objectDimTags=[(2, 50)], toolDimTags=[(2, 100), (2, 200), (2, 300), (2, 400)], tag=1000)
    # occ.rotate([(2, 1000)], 0, 0, 0, 0, 0, 1, math.pi/2)
    occ.synchronize()

    gmsh.model.addPhysicalGroup(2, [1000])
    # out = gmsh.model.occ.copy([(2, 1000)])
    # occ.translate(out, dx, 0, 0)
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

    # gmsh.initialize()


    setPeriodic(0)
    setPeriodic(1)
    gmsh.model.mesh.setSize(dimTags=gmsh.model.getEntities(0), size=lcar)
    gmsh.model.mesh.generate(dim=2)
    # gmsh.write('Hexg.msh')
    gmsh.write(mesh_filename + '.msh')
    os.system("gmsh -2 -o " + mesh_filename + ".msh -format msh22 " + mesh_filename + ".msh")
    os.system("dolfin-convert " + mesh_filename + ".msh " +  mesh_filename + ".xml")

    gmsh.write(mesh_filename+"-mesh.vtk")
    gmsh.finalize()

    mesh = meshio.read(mesh_filename+"-mesh.vtk")
    mesh.points = mesh.points[:, :2]
    meshio.write(mesh_filename+"-mesh.xdmf", mesh)
