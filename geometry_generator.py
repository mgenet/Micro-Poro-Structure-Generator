################################################################################
###                                                                          ###
### Created by Mahdi Manoochhertaybei, 2020-2023                             ###
###                                                                          ###
### École Polytechnique, Palaiseau, France                                   ###
###                                                                          ###
################################################################################

import geometry
import generate_seeds_random
import math
import os
import numpy

################################################################################

# os.system("gmsh Geometries/Hallow_square/Hex_rectangle2.geo -3 -format msh22 -o Geometries/Hallow_square/Hex_rectangle2.msh")
# os.system("dolfin-convert Geometries/Hallow_square/Hex_rectangle2.msh Geometries/Hallow_square/Hex_rectangle2.xml")

# choose the geometry
# geo = 'TKD'
# geo = 'hexagon_RVE1'
# geo = 'hexagon_RVE1_voronoi'
# geo = 'hexagon_RVE2'
# geo = 'Hallow_cube_inclusion_inside'
# geo = 'Hallow_cube_inclusion_corners'
# geo = 'Hallow_square_inclusion_inside'
# geo = 'Hallow_square_inclusion_corneres'
# geo = 'Hallow_square_inclusion_sides'
# geo = 'Hallow_array_inclusion_inside'
# geo = 'Hallow_array_inclusion_corneres'
# geo = 'Hallow_array_inclusion_sides'
geo = 'Hex_rectangle'
# geo = 'hexag_incl'
# geo = 'hexg'  # Hexagonal geometry has been created by substracting two hexagon from each others
# geo = 'Porous_media'
# geo= 'None'

Es = 1.0
nus = 0.2

if geo == 'TKD':
    l = 1
    h = 0.1

    fname = 'Geometries/TKD' + geo
    geometry.TKD(fname, l, h)

if geo == 'hexagon_RVE1':
    domain = 1
    row = 1
    DoI = 0.0
    thickness = 0.092
    shift_y = 0


    fname = 'Geometries/Hexagon/' + geo #+ '_voronoi'
    points = generate_seeds_random.semi_regular(DoI, row, domain)
    geometry.voronoi(fname, thickness, row, domain, shift_y, seeds_remove=True)

if geo == 'hexagon_RVE1_voronoi':

    domain = 1
    row = 3
    DoI = 0.4
    thickness = 0.17
    shift_y = 0

    fname = 'Geometries/Hexagon/' + geo #+ '_voronoi'
    # seeds.semi_regular(DoI, row, domain)
    # seeds.random(50)
    geometry.voronoi_tessellation(fname, thickness, row, domain, shift_y, seeds_remove=None)

if geo == 'hexagon_RVE2':
    domain = 1
    row = 1
    DoI = 0.0
    thickness = 0.1
    shift_y = -0.2

    fname = 'Geometries/Hexagon/' + geo
    generate_seeds_random.semi_regular(DoI, row, domain)
    geometry.voronoi(fname, thickness, row, domain, shift_y, seeds_remove=True)
    

if geo == 'Hallow_cube_inclusion_inside':

    phi = 0.7

    width = 1
    r = width/5
    # r = width * math.sqrt(phi/math.pi)
    
    shift_x = width/2
    shift_y = width/2
    shift_z = width/2

    fname = 'Geometries/Hallow_cube/' + geo
    geometry.Hallow_cube(fname, width, r, shift_x, shift_y, shift_z)


if geo == 'Hallow_cube_inclusion_corners':

    phi = 0.7

    width = 1
    r = width/5
    # r = width * math.sqrt(phi/math.pi)
    
    shift_x = 0
    shift_y = 0
    shift_z = 0

    fname = 'Geometries/Hallow_cube/' + geo
    geometry.Hallow_cube(fname, width, r, shift_x, shift_y, shift_z)


if geo == 'Hallow_square_inclusion_inside':

    phi = 0.1

    width = 0.25
    r = width/5
    # r = width * math.sqrt(phi/math.pi)
    
    shift_x = width/2
    shift_y = width/2

    fname = 'Geometries/Hallow_square/' + geo
    geometry.Hallow_Square(fname, width, r, shift_x, shift_y)

if geo == 'Hallow_square_inclusion_corneres':

    width = 0.25
    r = width/5
    # r = width * math.sqrt(phi/math.pi)

    shift_x = 0
    shift_y = 0

    fname = 'Geometries/Hallow_square/' + geo
    geometry.Hallow_Square(fname, width, r, shift_x, shift_y)

if geo == 'Hallow_square_inclusion_sides':

    width = 0.25
    r = width/5

    shift_x = width/2
    shift_y = 0

    fname = 'Geometries/Hallow_square/' + geo
    geometry.Hallow_Square(fname, width, r, shift_x, shift_y)




if geo == 'Hallow_array_inclusion_inside':

    width = 1
    shift_x = width/6
    shift_y = width/6

    fname = 'Geometries/Hallow_square/' + geo
    geometry.Hallow_array(fname, shift_x, shift_y, width)

if geo == 'Hallow_array_inclusion_corneres':

    width = 1
    shift_x = 0
    shift_y = 0

    fname = 'Geometries/Hallow_square/' + geo
    geometry.Hallow_array(fname, shift_x, shift_y, width)

if geo == 'Hallow_array_inclusion_sides':

    width = 1
    shift_x = width/6
    shift_y = 0

    fname = 'Geometries/Hallow_square/' + geo
    geometry.Hallow_array(fname, shift_x, shift_y, width)



if geo == 'Hex_rectangle':

    width = 1
    porosity = 0.5
    vol = math.sqrt(3)*width**2

    r = math.sqrt(porosity*vol/math.pi/2)
    # r = width/5

    # shift_x = 2*r
    # shift_y = 2*r

    shift_x = 0
    shift_y = 0

    fname = 'Geometries/Hallow_square/Hex_rectangle'
    geometry.Hex_rectangle(fname, width, r, shift_x, shift_y)


if geo == 'hexag_incl':

    width = 1
    angle = math.pi/3
    r = width/10
    # shift_x = 3*r
    # shift_y = 3*r
    shift_x = 0
    shift_y = 0

    fname = 'Geometries/Hexag_incl/hexag_incl'
    geometry.Hexag_incl(fname, width, angle, r, shift_x, shift_y)


    a = width     
    t = angle
    b = math.sin(t)*a
    c = a * math.cos(t)
    vol = a * b
    mesh_corners = numpy.array([[0, 0.],
                         [a, 0.],
                         [a+c, b],
                         [c, b]])

if geo == 'hexg':

    phi = 0.01
    R = 1

    fname = 'Geometries/Hexagon/Hexg'
    geometry.Hexg(fname, phi, R)



if geo == 'Porous_media':

    width = 0.25
    r = width/5

    shift_x = 0
    shift_y = 0

    fname = 'Geometries/Hallow_square/' + geo
    geometry.Porous_media(fname, width, r, shift_x, shift_y)

# if fname == 'Geometries/Hexag_incl/hexag_incl':
#     homogenized_parameters.homogenization(fname, Es, nus, mesh_corners, vol)
# else:
#     homogenized_parameters.homogenization(fname, Es, nus)
# homogenized_parameters.homogenization('Geometries/Hallow_square/Hex_rectangle2', Es, nus)


