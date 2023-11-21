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

def generate_mesh_2D_parallelogram_w_circular_inclusions(
        mesh_filename,
        width,
        angle,
        r,
        lcar,
        shift_x = 0.,
        shift_y = 0.):

    ############# Parameters ###################################################

    a = width
    R = r
    t = angle
    b = math.sin(t)*a
    c = a * math.cos(t)

    xmin = 0; xmax = c
    ymin = 0; ymax = b
    zmin = 0; zmax = 0
    e = 1e-6

    ############# Generator ####################################################

    gmsh.initialize()
    gmsh.clear()
    model = gmsh.model
    occ = model.occ

    occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(occ.addPoint(0, 0, 0, lcar), occ.addPoint(a, 0, 0, lcar)),\
                                           occ.addLine(occ.addPoint(a, 0, 0, lcar), occ.addPoint(a+c, b, 0, lcar)),\
                                           occ.addLine(occ.addPoint(a+c, b, 0, lcar), occ.addPoint(c, b, 0, lcar)),\
                                           occ.addLine(occ.addPoint(c, b, 0, lcar), occ.addPoint(0, 0, 0, lcar))])], 1)

    occ.addDisk(xc=0+shift_x, yc=0+shift_y, zc=0, rx=R, ry=R, tag=2)
    occ.addDisk(xc=a+shift_x, yc=0+shift_y, zc=0, rx=R, ry=R, tag=3)
    occ.addDisk(xc=a+c+shift_x, yc=b+shift_y, zc=0, rx=R, ry=R, tag=4)
    occ.addDisk(xc=c+shift_x, yc=b+shift_y, zc=0, rx=R, ry=R, tag=5)
    occ.cut(objectDimTags=[(2, 1)], toolDimTags=[(2, 2), (2, 3), (2, 4), (2, 5)], tag=6)
    occ.synchronize()
    gmsh.model.addPhysicalGroup(2, [6])

    gen.setPeriodic(coord=0, xmin=xmin, ymin=ymin, zmin=zmin, xmax=xmax, ymax=ymax, zmax=zmax, e=e)
    gen.setPeriodic(coord=1, xmin=xmin, ymin=ymin, zmin=zmin, xmax=xmax, ymax=ymax, zmax=zmax, e=e)
    gmsh.model.mesh.setSize(dimTags=gmsh.model.getEntities(0), size=lcar)
    gmsh.model.mesh.generate(dim=2)

    gmsh.write(mesh_filename + '.msh')

    gen.convert_msh_to_xml(mesh_filename)

    gmsh.write(mesh_filename+".vtk")

    gen.convert_vtk_to_xdmf(mesh_filename, dim=2)

    gmsh.finalize()
