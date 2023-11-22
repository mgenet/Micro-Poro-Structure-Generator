################################################################################
###                                                                          ###
### Created by Mahdi Manoochehrtayebi, 2020-2023                             ###
###                                                                          ###
### École Polytechnique, Palaiseau, France                                   ###
###                                                                          ###
################################################################################

import math
import gmsh

import micro_poro_structure_generator as gen

################################################################################

def generate_mesh_3D_cube_w_spherical_inclusions(
        mesh_filename,
        width,
        r,
        lcar,
        shift_x = 0.,
        shift_y = 0.,
        shift_z = 0.):

    ########### Initialization #################################################
    
    xmin = 0.
    ymin = 0.
    zmin = 0.
    xmax = width
    ymax = width
    zmax = width

    r0 = r
    e = 1e-6

    ############# Generator ####################################################

    gmsh.initialize()

    box_tag = 1
    hole_tag1 = 2
    hole_tag2 = 3
    hole_tag3 = 4
    hole_tag4 = 5
    hole_tag5 = 6
    hole_tag6 = 7
    hole_tag7 = 8
    hole_tag8 = 9
    rve_tag = 10

    gmsh.model.occ.addBox(x=xmin, y=ymin, z=zmin, dx=xmax-xmin, dy=ymax-ymin, dz=zmax-zmin, tag=box_tag)
    gmsh.model.occ.addSphere(xc=xmin+shift_x, yc=ymin+shift_y, zc=zmin+shift_z, radius=r0, tag=hole_tag1)
    gmsh.model.occ.addSphere(xc=xmax+shift_x, yc=ymin+shift_y, zc=zmin+shift_z, radius=r0, tag=hole_tag2)
    gmsh.model.occ.addSphere(xc=xmax+shift_x, yc=ymax+shift_y, zc=zmin+shift_z, radius=r0, tag=hole_tag3)
    gmsh.model.occ.addSphere(xc=xmin+shift_x, yc=ymax+shift_y, zc=zmin+shift_z, radius=r0, tag=hole_tag4)
    gmsh.model.occ.addSphere(xc=xmin+shift_x, yc=ymin+shift_y, zc=zmax+shift_z, radius=r0, tag=hole_tag5)
    gmsh.model.occ.addSphere(xc=xmax+shift_x, yc=ymin+shift_y, zc=zmax+shift_z, radius=r0, tag=hole_tag6)
    gmsh.model.occ.addSphere(xc=xmax+shift_x, yc=ymax+shift_y, zc=zmax+shift_z, radius=r0, tag=hole_tag7)
    gmsh.model.occ.addSphere(xc=xmin+shift_x, yc=ymax+shift_y, zc=zmax+shift_z, radius=r0, tag=hole_tag8)
    gmsh.model.occ.cut(objectDimTags=[(3, box_tag)], toolDimTags=[(3, hole_tag1), (3, hole_tag2), (3, hole_tag3), (3, hole_tag4), (3, hole_tag5), (3, hole_tag6), (3, hole_tag7), (3, hole_tag8)], tag=rve_tag)
    gmsh.model.occ.synchronize()
    gmsh.model.addPhysicalGroup(dim=3, tags=[rve_tag])
    gen.setPeriodic(coord=0, xmin=xmin, ymin=ymin, zmin=zmin, xmax=xmax, ymax=ymax, zmax=zmax, e=e)
    gen.setPeriodic(coord=1, xmin=xmin, ymin=ymin, zmin=zmin, xmax=xmax, ymax=ymax, zmax=zmax, e=e)
    gen.setPeriodic(coord=2, xmin=xmin, ymin=ymin, zmin=zmin, xmax=xmax, ymax=ymax, zmax=zmax, e=e)

    gmsh.model.mesh.setSize(dimTags=gmsh.model.getEntities(0), size=lcar)
    gmsh.model.mesh.generate(dim=3)

    gmsh.write(mesh_filename + '.msh')

    gen.convert_msh_to_xml(mesh_filename)

    gmsh.write(mesh_filename+".vtk")

    gen.convert_vtk_to_xdmf(mesh_filename, dim=3)

    gmsh.finalize()

    vol = width**3
    phi = (4/3*math.pi*r**3)/vol
    print("porosity:" +str(phi))
    return phi
