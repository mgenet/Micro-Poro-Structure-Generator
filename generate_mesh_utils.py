################################################################################
###                                                                          ###
### Created by Mahdi Manoochehrtayebi, 2020-2023                             ###
###                                                                          ###
### École Polytechnique, Palaiseau, France                                   ###
###                                                                          ###
################################################################################

import gmsh
import os

################################################################################

def setPeriodic(coord, xmin, ymin, zmin, xmax, ymax, zmax, e):
    # From https://gitlab.onelab.info/gmsh/gmsh/-/issues/744

    smin = gmsh.model.getEntitiesInBoundingBox(
        xmin - e,
        ymin - e,
        zmin - e,
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
            bb2[0] -= dx; bb2[1] -= dy; bb2[2] -= dz
            bb2[3] -= dx; bb2[4] -= dy; bb2[5] -= dz
            if ((abs(bb2[0] - bb[0]) < e) and (abs(bb2[1] - bb[1]) < e) and
                (abs(bb2[2] - bb[2]) < e) and (abs(bb2[3] - bb[3]) < e) and
                (abs(bb2[4] - bb[4]) < e) and (abs(bb2[5] - bb[5]) < e)):
                gmsh.model.mesh.setPeriodic(2, [j[1]], [i[1]], [1, 0, 0, dx,\
                                                                0, 1, 0, dy,\
                                                                0, 0, 1, dz,\
                                                                0, 0, 0, 1 ])

################################################################################

def compute_porosity_2D_using_fenics(mesh_filename):

    import dolfin

    mesh = dolfin.Mesh()
    dolfin.XDMFFile(mesh_filename+".xdmf").read(mesh)
    coord = mesh.coordinates()
    xmax = max(coord[:,0]); xmin = min(coord[:,0])
    ymax = max(coord[:,1]); ymin = min(coord[:,1])
    V0 = (xmax - xmin)*(ymax - ymin)
    dV = dolfin.Measure("dx", domain=mesh)
    Vs0 = dolfin.assemble(dolfin.Constant(1.) * dV)

    Phis0 = Vs0/V0
    Phif0 = 1. - Phis0

    print ("Phis0=" +str(Phis0))
    print ("Phif0=" +str(Phif0))

    return Phif0

################################################################################

def convert_msh_to_xml(mesh_filename):

    os.system("gmsh -2 -o " + mesh_filename + ".msh -format msh22 " + mesh_filename + ".msh")
    os.system("dolfin-convert " + mesh_filename + ".msh " + mesh_filename + ".xml")

################################################################################

def convert_vtk_to_xdmf(mesh_filename, dim=3):

    import meshio

    mesh = meshio.read(mesh_filename+".vtk")
    if (dim == 2):
        mesh.points = mesh.points[:, :2]
    meshio.write(mesh_filename+".xdmf", mesh)
