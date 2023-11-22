################################################################################
###                                                                          ###
### Created by Mahdi Manoochehrtayebi, 2020-2023                             ###
###                                                                          ###
### École Polytechnique, Palaiseau, France                                   ###
###                                                                          ###
################################################################################

import gmsh
import math

import micro_poro_structure_generator as gen

################################################################################

def generate_mesh_2D_rectangle_w_hexagonal_inclusions(
        mesh_filename,
        phi,
        R,
        lcar):

    ############# Parameters ###################################################

    h = R*(1 - phi)
    L = R + h
    lcar = h/10

    ############# Generator ####################################################

    gmsh.initialize()
    gmsh.clear()
    model = gmsh.model
    occ = model.occ

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

    hexagon_generator(0, 0, 10)
    hexagon_generator((2*R + h)*math.sqrt(3)/2, 0, 20)
    hexagon_generator((2*R + h)*math.sqrt(3)/4, (2*R + h)*3/4, 30)
    hexagon_generator((2*R + h)*math.sqrt(3)/4, -(2*R + h)*3/4, 40)
    occ.fuse([(2,12)], [(2,22), (2, 32), (2, 42)], 50)
    shift = 0.
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

    gen.setPeriodic(coord=0, xmin=xmin, ymin=ymin, zmin=zmin, xmax=xmax, ymax=ymax, zmax=zmax, e=e)
    gen.setPeriodic(coord=1, xmin=xmin, ymin=ymin, zmin=zmin, xmax=xmax, ymax=ymax, zmax=zmax, e=e)
    gmsh.model.mesh.setSize(dimTags=gmsh.model.getEntities(0), size=lcar)
    gmsh.model.mesh.generate(dim=2)

    gmsh.write(mesh_filename + '.msh')

    gen.convert_msh_to_xml(mesh_filename)

    gmsh.write(mesh_filename+".vtk")

    gen.convert_vtk_to_xdmf(mesh_filename, dim=2)

    gmsh.finalize()
