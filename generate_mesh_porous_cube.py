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

def generate_mesh_porous_cube(mesh_filename, width, r, shift_x, shift_y, shift_z):


    ########### Initialization #################################################
    
    xmin = 0.
    ymin = 0.
    zmin = 0.
    xmax = width
    ymax = width
    zmax = width

    r0 = r
    l = width/10
    e = 1e-6


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
    setPeriodic(0)
    setPeriodic(1)
    setPeriodic(2)

    gmsh.model.mesh.setSize(dimTags=gmsh.model.getEntities(0), size=l)
    gmsh.model.mesh.generate(dim=3)

    # gmsh.model.mesh.generate()
    gmsh.write(mesh_filename + '.msh')
    os.system("gmsh -2 -o " + mesh_filename + ".msh -format msh22 " + mesh_filename + ".msh")
    os.system("dolfin-convert " + mesh_filename + ".msh " +  mesh_filename + ".xml")

    gmsh.write(mesh_filename+"-mesh.vtk")
    gmsh.finalize()

    mesh = meshio.read(mesh_filename+"-mesh.vtk")
    mesh.points = mesh.points[:, :3]
    meshio.write(mesh_filename+"-mesh.xdmf", mesh)

    mesh = dolfin.Mesh()
    dolfin.XDMFFile(mesh_filename+"-mesh.xdmf").read(mesh)
    
    vol = width**3
    phi = (4/3*math.pi*r**3)/vol
    print("porosity:" +str(phi))
