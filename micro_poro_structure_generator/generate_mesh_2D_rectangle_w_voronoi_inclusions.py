################################################################################
###                                                                          ###
### Created by Mahdi Manoochehrtayebi, 2020-2023                             ###
###                                                                          ###
### École Polytechnique, Palaiseau, France                                   ###
###                                                                          ###
################################################################################

import gmsh
import math
import numpy
import os
import pickle
import scipy
import scipy.spatial

import micro_poro_structure_generator as gen

################################################################################

def generate_mesh_2D_rectangle_w_voronoi_inclusions(
        mesh_filename,
        seeds_filename,
        h,                   # wall thickness
        lcar,                # element size
        domain_x,
        domain_y,
        shift_y = 0.,
        remove_seeds = True):

    ############# Functions ####################################################

    def line_btw_points(P, Q):
        # Returns the line parameters between two points
        
        a = Q[1] - P[1]
        b = Q[0] - P[0]
        c = b * P[1] - a * P[0]
        return a, b, c

    def line_2_points(P, Q):
        # Returns the line parameters between two points
        
        m = (Q[1] - P[1])/(Q[0] - P[0])
        c = P[1] - m * P[0]
        return m, c

    def perpendicular_line(P, Q, a, b, c):
        # Returns the perpendicular line to a line passing two points
        
        mid_point = [(P[0] + Q[0])/2, (P[1] + Q[1])/2]
        c = b * (mid_point[0]) + a * (mid_point[1])
        temp = a
        a = -b
        b = temp
        return a, b, c

    def lines_intersect(m1, c1, m2, c2):
        # Returns the intersection point of two lines
        
        x = (c2 - c1)/(m1 - m2)
        y = (m2*c1 - m1*c2)/(m2 - m1)
        return [x, y]

    def lines_intersection(a1, b1, c1, a2, b2, c2):
        # Returns the intersection point of two lines
        
        determinant = a1 * b2 - a2 * b1
        if (determinant == 0):
            # The lines are parallel. This is simplified
            # by returning a pair of (10.0)**19
            return [(10.0)**19, (10.0)**19]
        else:
            x = (b1 * c2 - b2 * c1)/determinant
            y = (a1 * c2 - a2 * c1)/determinant
            return [x, y]

    def vertex_generator(P, Q, S):

        a1, b1, c1 = line_btw_points(P, Q)
        a2, b2, c2 = line_btw_points(Q, S)
        a1, b1, c1 = perpendicular_line(P, Q, a1,b1,c1)
        a2, b2, c2 = perpendicular_line(Q, S, a2,b2,c2)
        x, y = lines_intersection(a1,b1,c1,a2,b2,c2)
        return [x, y]

    def find_neighbor_triangles(triangle_num, triangle, triangles_cor):

        neighbors = []
        neighbors.append(triangle_num)
        for i in range(len(triangles_cor)):
            relation = 0
            for j in range(3):
                for k in range(3):
                    if triangle[k] == triangles_cor[i][j]:
                        relation += 1
            if relation == 2:
                neighbors.append(i)

        return neighbors

    def one_points_intersection(neighbor_vertices, vertices, lines):

        lines.append([vertices[neighbor_vertices[0]][1]])
        return lines

    def two_points_intersection(neighbor_vertices, vertices, lines):

        mid_point = vertices[neighbor_vertices[0]][1]
        side_point_1 = vertices[neighbor_vertices[1]][1]
        side_point_2 = vertices[neighbor_vertices[2]][1]

        mid_point_num = neighbor_vertices[0]
        side_point_1_num = neighbor_vertices[1]
        side_point_2_num = neighbor_vertices[2]

        m1, c1 = line_2_points(mid_point,side_point_1)
        beta = math.atan(m1)
        sin = math.sin(beta)
        cos = math.cos(beta)

        if side_point_1[0] < mid_point[0]:
            point_up_1 = [mid_point[0] - h/2 * sin, mid_point[1] + h/2 * cos]
            c1_up = point_up_1[1] - m1 * point_up_1[0]
            point_down_1 = [mid_point[0] + h/2 * sin, mid_point[1] - h/2 * cos]
            c1_down = point_down_1[1] - m1 * point_down_1[0]

        else:
            point_down_1 = [mid_point[0] - h/2 * sin, mid_point[1] + h/2 * cos]
            c1_down = point_down_1[1] - m1 * point_down_1[0]
            point_up_1 = [mid_point[0] + h/2 * sin, mid_point[1] - h/2 * cos]
            c1_up = point_up_1[1] - m1 * point_up_1[0]

        m2, c2 = line_2_points(mid_point,side_point_2)
        beta = math.atan(m2)
        sin = math.sin(beta)
        cos = math.cos(beta)

        if side_point_2[0] < mid_point[0]:
            point_up_2 = [mid_point[0] - h/2 * sin, mid_point[1] + h/2 * cos]
            c2_up = point_up_2[1] - m2 * point_up_2[0]
            point_down_2 = [mid_point[0] + h/2 * sin, mid_point[1] - h/2 * cos]
            c2_down = point_down_2[1] - m2 * point_down_2[0]

        else:
            point_down_2 = [mid_point[0] - h/2 * sin, mid_point[1] + h/2 * cos]
            c2_down = point_down_2[1] - m2 * point_down_2[0]
            point_up_2 = [mid_point[0] + h/2 * sin, mid_point[1] - h/2 * cos]
            c2_up = point_up_2[1] - m2 * point_up_2[0]

        mid_point_up = lines_intersect(m1, c1_up, m2, c2_down)
        mid_point_down = lines_intersect(m1, c1_down, m2, c2_up)

        lines.append([mid_point_num, mid_point_up, mid_point_down])

        return lines

    def getAngle(a, b, c):

        ang = math.degrees(math.atan2(c[1]-b[1], c[0]-b[0]) - math.atan2(a[1]-b[1], a[0]-b[0]))
        return ang%360

    def three_points_intersection(neighbor_vertices, vertices, lines):
        # This function gets a triangle center seed which has three neighbors and returns a small mesh triangle around the triangle centers

        mid_point = vertices[neighbor_vertices[0]][1]
        side_point_1 = vertices[neighbor_vertices[1]][1]
        side_point_2 = vertices[neighbor_vertices[2]][1]
        side_point_3 = vertices[neighbor_vertices[3]][1]
        side_points = [side_point_1, side_point_2, side_point_3]

        mid_point_num = neighbor_vertices[0]
        side_point_1_num = neighbor_vertices[1]
        side_point_2_num = neighbor_vertices[2]
        side_point_3_num = neighbor_vertices[3]

        m = numpy.zeros(3)
        c = numpy.zeros(3)
        c_up = numpy.zeros(3)
        c_down = numpy.zeros(3)
        for i in range(len(side_points)):
            m[i], c[i] = line_2_points(mid_point, side_points[i])
            beta = math.atan(m[i])
            sin = math.sin(beta)
            cos = math.cos(beta)

            if side_points[i][0] < mid_point[0]:
                point_up = [mid_point[0] - h/2 * sin, mid_point[1] + h/2 * cos]
                c_up[i] = point_up[1] - m[i] * point_up[0]
                point_down = [mid_point[0] + h/2 * sin, mid_point[1] - h/2 * cos]
                c_down[i] = point_down[1] - m[i] * point_down[0]

            else:
                point_down = [mid_point[0] - h/2 * sin, mid_point[1] + h/2 * cos]
                c_down[i] = point_down[1] - m[i] * point_down[0]
                point_up = [mid_point[0] + h/2 * sin, mid_point[1] - h/2 * cos]
                c_up[i] = point_up[1] - m[i] * point_up[0]

        gamma2 = getAngle(side_points[0], mid_point, side_points[1])
        gamma3 = getAngle(side_points[0], mid_point, side_points[2])

        # If clockwise:
        if gamma2 > gamma3:
            mid_point_12 = lines_intersect(m[0], c_up[0], m[1], c_down[1])
            mid_point_23 = lines_intersect(m[1], c_up[1], m[2], c_down[2])
            mid_point_31 = lines_intersect(m[2], c_up[2], m[0], c_down[0])

        # If counter-clockwise:
        else:
            mid_point_12 = lines_intersect(m[0], c_down[0], m[1], c_up[1])
            mid_point_23 = lines_intersect(m[1], c_down[1], m[2], c_up[2])
            mid_point_31 = lines_intersect(m[2], c_down[2], m[0], c_up[0])

        lines.append([mid_point_num, [side_point_1_num , side_point_2_num, mid_point_12], [side_point_2_num, side_point_3_num, mid_point_23], [side_point_3_num, side_point_1_num, mid_point_31]])

        return lines

    ############# Generator ####################################################

    gmsh.initialize()
    gmsh.clear()
    model = gmsh.model
    occ = model.occ

    open_file = open(seeds_filename, "rb")
    points = pickle.load(open_file)
    open_file.close()

    if (remove_seeds is True):
        os.remove(seeds_filename)

    nb_points = len(points)

    # vor = scipy.spatial.Voronoi(points)
    # fig = scipy.spatial.voronoi_plot_2d(vor)
    # plt.savefig('Voronoi.jpg')

    period_neighbor_points_1 = numpy.zeros((nb_points,2))
    period_neighbor_points_2 = numpy.zeros((nb_points,2))
    period_neighbor_points_3 = numpy.zeros((nb_points,2))
    period_neighbor_points_4 = numpy.zeros((nb_points,2))
    period_neighbor_points_5 = numpy.zeros((nb_points,2))
    period_neighbor_points_6 = numpy.zeros((nb_points,2))
    period_neighbor_points_7 = numpy.zeros((nb_points,2))
    period_neighbor_points_8 = numpy.zeros((nb_points,2))

    for i in range(len(points)):
        period_neighbor_points_1[i][0] = points[i][0] - domain_x
        period_neighbor_points_1[i][1] = points[i][1] + domain_y

        period_neighbor_points_2[i][0] = points[i][0]
        period_neighbor_points_2[i][1] = points[i][1] + domain_y

        period_neighbor_points_3[i][0] = points[i][0] + domain_x
        period_neighbor_points_3[i][1] = points[i][1] + domain_y

        period_neighbor_points_4[i][0] = points[i][0] - domain_x
        period_neighbor_points_4[i][1] = points[i][1]

        period_neighbor_points_5[i][0] = points[i][0] + domain_x
        period_neighbor_points_5[i][1] = points[i][1]

        period_neighbor_points_6[i][0] = points[i][0] - domain_x
        period_neighbor_points_6[i][1] = points[i][1] - domain_y

        period_neighbor_points_7[i][0] = points[i][0]
        period_neighbor_points_7[i][1] = points[i][1] - domain_y

        period_neighbor_points_8[i][0] = points[i][0] + domain_x
        period_neighbor_points_8[i][1] = points[i][1] - domain_y
    
    points = numpy.concatenate((points, period_neighbor_points_1, period_neighbor_points_2, period_neighbor_points_3, period_neighbor_points_4, period_neighbor_points_5, period_neighbor_points_6, period_neighbor_points_7, period_neighbor_points_8))

    tri = scipy.spatial.Delaunay(points)

    # plt.plot(points[:,0], points[:,1], 'o')
    # plt.triplot(points[:,0], points[:,1], tri.simplices)

    # Put each triangle corners in one group
    triangles_cor = tri.simplices
    number_of_triangles = len(triangles_cor)

    # Find vertices for each triangle
    # Put each triangle corners and its vertices in one group and make a list, please
    vertices = []
    for i in range(number_of_triangles):
        vertices.append([triangles_cor[i], vertex_generator(points[triangles_cor[i, 0]], points[triangles_cor[i, 1]], points[triangles_cor[i, 2]])])

    # Find the neighbor triangles
    neighbor_triangles_by_number = []
    for i in range(number_of_triangles):
        neighbor_triangles_by_number.append(find_neighbor_triangles(i, triangles_cor[i], triangles_cor))

    # Put neighbor triangles vertices in one group which shows the ridge of Voronoi tessellation
    lines = []      # lines contains each triangle (vertices) number with the side points
    for i in range (len(neighbor_triangles_by_number)):
        if len(neighbor_triangles_by_number[i]) == 2:
            lines = one_points_intersection(neighbor_triangles_by_number[i], vertices, lines)
        if len(neighbor_triangles_by_number[i]) == 3:
            lines = two_points_intersection(neighbor_triangles_by_number[i], vertices, lines)
        if len(neighbor_triangles_by_number[i]) == 4:
            lines = three_points_intersection(neighbor_triangles_by_number[i], vertices, lines)

    surfaces = []
    l = 0
    for i in range(len(neighbor_triangles_by_number)):
        if len(neighbor_triangles_by_number[i]) == 4:
            point = vertices[i][1]
            for j in range(len(neighbor_triangles_by_number[i])-1):
                neighbor_num = neighbor_triangles_by_number[i][j+1]
                if neighbor_num > i:
                    neighbor = vertices[neighbor_num][1]
                    m, c = line_2_points(point, neighbor)
                    beta = math.atan(m)
                    sin = math.sin(beta)
                    cos = math.cos(beta)

                    point_up = [point[0] - h/2 * sin, point[1] + h/2 * cos]
                    point_down = [point[0] + h/2 * sin, point[1] - h/2 * cos]

                    neighbor_up = [neighbor[0] - h/2 * sin, neighbor[1] + h/2 * cos]
                    neighbor_down = [neighbor[0] + h/2 * sin, neighbor[1] - h/2 * cos]

                    p_up = occ.addPoint(point_up[0], point_up[1], 0, lcar)
                    p_down = occ.addPoint(point_down[0], point_down[1], 0, lcar)
                    n_up = occ.addPoint(neighbor_up[0], neighbor_up[1], 0, lcar)
                    n_down = occ.addPoint(neighbor_down[0], neighbor_down[1], 0, lcar)

                    l += 1
                    if point[0] > neighbor[0]:
                        occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(p_up, p_down), occ.addLine(p_down, n_down), occ.addLine(n_down, n_up), occ.addLine(n_up, p_up)])], l)

                    else:
                        occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(p_down, p_up), occ.addLine(p_up, n_up), occ.addLine(n_up, n_down), occ.addLine(n_down, p_down)])], l)

                    surfaces.append([i, neighbor_num, l])

    occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(occ.addPoint(0, 0+shift_y, 0, lcar), occ.addPoint(domain_x, 0+shift_y, 0, lcar)),\
                                           occ.addLine(occ.addPoint(domain_x, 0+shift_y, 0, lcar), occ.addPoint(domain_x, domain_y+shift_y, 0, lcar)),\
                                           occ.addLine(occ.addPoint(domain_x, domain_y+shift_y, 0, lcar), occ.addPoint(0, domain_y+shift_y, 0, lcar)),\
                                           occ.addLine(occ.addPoint(0, domain_y+shift_y, 0, lcar), occ.addPoint(0, 0+shift_y, 0, lcar))])], 2000000)

    base = 'occ.fuse([(2,1)], [(2,2)'
    for i in range(l-2):
        base = base + ',(2,' + str(i+3) + ')'
    command = base + '], 3000000)'
    exec(command)

    frame = occ.cut([(2, 2000000)], [(2, 3000000)])

    occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(occ.addPoint(0, 0+shift_y, 0, lcar), occ.addPoint(domain_x, 0+shift_y, 0, lcar)),\
                                           occ.addLine(occ.addPoint(domain_x, 0+shift_y, 0, lcar), occ.addPoint(domain_x, domain_y+shift_y, 0, lcar)),\
                                           occ.addLine(occ.addPoint(domain_x, domain_y+shift_y, 0, lcar), occ.addPoint(0, domain_y+shift_y, 0, lcar)),\
                                           occ.addLine(occ.addPoint(0, domain_y+shift_y, 0, lcar), occ.addPoint(0, 0+shift_y, 0, lcar))])], 2000001)

    base = 'occ.cut([(2, 2000001)], [(2, 1)'
    for i in range(len(frame[0])-1):
        base = base + ',(2,' + str(i+2) + ')'
    command = base + '], 4000000)'
    exec(command)

    gmsh.option.setNumber("Mesh.MeshSizeMax", lcar)
    
    occ.synchronize()
    gmsh.model.addPhysicalGroup(2, [4000000])
    
    zmin = 0
    zmax = 0
    ymin = 0+shift_y
    ymax = domain_y+shift_y
    xmin = 0
    xmax = domain_x
    e = 1e-6
    gen.setPeriodic(coord=0, xmin=xmin, ymin=ymin, zmin=zmin, xmax=xmax, ymax=ymax, zmax=zmax, e=e)
    gen.setPeriodic(coord=1, xmin=xmin, ymin=ymin, zmin=zmin, xmax=xmax, ymax=ymax, zmax=zmax, e=e)

    gmsh.model.mesh.generate()

    gmsh.write(mesh_filename+".msh")
    gmsh.write(mesh_filename+".vtk")
    gmsh.finalize()

    gen.convert_msh_to_xml(mesh_filename)

    gen.convert_vtk_to_xdmf(mesh_filename, dim=2)

    Phif0 = gen.compute_porosity_2D_using_fenics(mesh_filename)
    return Phif0
