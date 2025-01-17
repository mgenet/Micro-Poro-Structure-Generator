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

def setPeriodic(dim, coord, xmin, ymin, zmin, xmax, ymax, zmax, e=1e-6):
    # From https://gitlab.onelab.info/gmsh/gmsh/-/issues/744

    dx = (xmax - xmin) if (coord == 0) else 0.
    dy = (ymax - ymin) if (coord == 1) else 0.
    dz = (zmax - zmin) if (coord == 2) else 0.
    d = max(dx, dy, dz)
    e *= d

    smin = gmsh.model.getEntitiesInBoundingBox(
        xmin      - e, ymin      - e, zmin      - e,
        xmax - dx + e, ymax - dy + e, zmax - dz + e,
        dim-1)
    # print ("smin:",smin)
    for i in smin:
        bb = gmsh.model.getBoundingBox(*i)
        bbe = [bb[0] + dx, bb[1] + dy, bb[2] + dz,
               bb[3] + dx, bb[4] + dy, bb[5] + dz]
        smax = gmsh.model.getEntitiesInBoundingBox(
            bbe[0] - e, bbe[1] - e, bbe[2] - e,
            bbe[3] + e, bbe[4] + e, bbe[5] + e,
            dim-1)
        # print ("smax:",smax)
        for j in smax:
            bb2 = gmsh.model.getBoundingBox(*j)
            bb2e = [bb2[0] - dx, bb2[1] - dy, bb2[2] - dz,
                    bb2[3] - dx, bb2[4] - dy, bb2[5] - dz]
            if (numpy.linalg.norm(numpy.asarray(bb2e) - numpy.asarray(bb)) < e):
                gmsh.model.mesh.setPeriodic(
                    dim-1,
                    [j[1]], [i[1]],
                    [1, 0, 0, dx,\
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
