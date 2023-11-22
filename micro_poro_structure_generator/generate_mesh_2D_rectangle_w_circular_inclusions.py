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

def generate_mesh_2D_rectangle_w_circular_inclusions(
        mesh_filename,
        width,
        r,
        lcar,
        shift_x = 0.,
        shift_y = 0.):
    
    ############# Parameters ###################################################

    length = width * math.sqrt(3)
    xmin = 0.
    ymin = 0.
    zmin = 0.
    xmax = width
    ymax = length
    zmax = 0
    x0 = width/2
    y0 = length/2
    z0 = 0
    r0 = r
    e = 1e-6

    ############# Generator ####################################################

    gmsh.initialize()
    gmsh.clear()

    box_tag = 1
    hole_tag = 2
    hole_tag2 = 3
    hole_tag3 = 4
    hole_tag4 = 5
    hole_tag5 = 6
    hole_tag6 = 7
    rve_tag = 8

    gmsh.model.occ.addRectangle(x=xmin+shift_x, y=ymin+shift_y, z=0, dx=xmax-xmin, dy=ymax-ymin, tag=box_tag)
    gmsh.model.occ.addDisk(xc=x0, yc=y0, zc=0, rx=r0, ry=r0, tag=hole_tag)
    gmsh.model.occ.addDisk(xc=xmin, yc=ymin, zc=0, rx=r0, ry=r0, tag=hole_tag2)
    gmsh.model.occ.addDisk(xc=xmax, yc=ymin, zc=0, rx=r0, ry=r0, tag=hole_tag3)
    gmsh.model.occ.addDisk(xc=xmax, yc=ymax, zc=0, rx=r0, ry=r0, tag=hole_tag4)
    gmsh.model.occ.addDisk(xc=xmin, yc=ymax, zc=0, rx=r0, ry=r0, tag=hole_tag5)
    gmsh.model.occ.addDisk(xc=x0 + xmax - xmin, yc=y0, zc=0, rx=r0, ry=r0, tag=hole_tag6)
    gmsh.model.occ.cut(objectDimTags=[(2, box_tag)], toolDimTags=[(2, hole_tag), (2, hole_tag2), (2, hole_tag3), (2, hole_tag4), (2, hole_tag5), (2, hole_tag6)], tag=rve_tag)
    gmsh.model.occ.synchronize()
    gmsh.model.addPhysicalGroup(dim=2, tags=[rve_tag])
    
    gen.setPeriodic(coord=0, xmin=xmin, ymin=ymin, zmin=zmin, xmax=xmax, ymax=ymax, zmax=zmax, e=e)
    gen.setPeriodic(coord=1, xmin=xmin, ymin=ymin, zmin=zmin, xmax=xmax, ymax=ymax, zmax=zmax, e=e)
    gmsh.model.mesh.setSize(dimTags=gmsh.model.getEntities(0), size=lcar)
    gmsh.model.mesh.generate(dim=2)

    gmsh.write(mesh_filename + '.msh')

    gen.convert_msh_to_xml(mesh_filename)

    gmsh.write(mesh_filename+".vtk")

    gen.convert_vtk_to_xdmf(mesh_filename, dim=2)

    gmsh.finalize()

    Phif0 = gen.compute_porosity_2D_using_fenics(mesh_filename)
    return Phif0
