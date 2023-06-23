################################################################################
###                                                                          ###
### Created by Mahdi Manoochhertaybei, 2020-2023                             ###
###                                                                          ###
### École Polytechnique, Palaiseau, France                                   ###
###                                                                          ###
################################################################################


def voronoi(fname, h, row, domain_y, shift_y, seeds_remove = True):

    import math
    import meshio
    import gmsh
    import os
    import pickle
    import dolfin
    import numpy as np
    from scipy.spatial import Delaunay
    import matplotlib.pyplot as plt
    from scipy.spatial import Voronoi, voronoi_plot_2d


    ########## Initialization#######################################################


    gmsh.initialize()
    gmsh.clear()
    model = gmsh.model
    occ = model.occ

    ############# Parameters #######################################################

    lcar = h/5
    l = 0
    ############# Funtions #########################################################

    def line_btw_points(P, Q):
    #Retunrs the line parameter between two points
        # by = ax + c
        a = Q[1] - P[1]
        b = Q[0] - P[0]
        c = b * P[1] - a * P[0]
        # beta = math.atan(a/b)
        return a, b, c

    def line_2_points(P, Q):
    #Retunrs the line parameter between two points
        m = (Q[1] - P[1])/(Q[0] - P[0])
        c = P[1] - m * P[0]
        return m, c

    def perpendicular_line(P, Q, a, b, c):
    #Returns the perpendicular line to a line passing two points
        mid_point = [(P[0] + Q[0])/2, (P[1] + Q[1])/2]

        # ay = -bx + d
        c = b * (mid_point[0]) + a * (mid_point[1])
        temp = a
        a = -b
        b = temp
        return a, b, c

    def lines_Intersect(m1, c1, m2, c2):
    # Returns the intersection point of two lines
        x = (c2 - c1)/(m1 - m2)
        y = (m2*c1 - m1*c2)/(m2 - m1)
        return [x, y]

    def lines_Intersection(a1, b1, c1, a2, b2, c2):
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


    def vertice_generator(P, Q, S):
        a1, b1, c1 = line_btw_points(P, Q)
        a2, b2, c2 = line_btw_points(Q, S)
        a1, b1, c1 = perpendicular_line(P, Q, a1,b1,c1)
        a2, b2, c2 = perpendicular_line(Q, S, a2,b2,c2)
        x, y = lines_Intersection(a1,b1,c1,a2,b2,c2)
        return [x, y]

    def find_neighbour_triangles(triangle_num, triangle, triangles_cor):
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




        mid_point_up = lines_Intersect(m1, c1_up, m2, c2_down)
        mid_point_down = lines_Intersect(m1, c1_down, m2, c2_up)

        lines.append([mid_point_num, mid_point_up, mid_point_down])
        return lines

    def getAngle(a, b, c):
        ang = math.degrees(math.atan2(c[1]-b[1], c[0]-b[0]) - math.atan2(a[1]-b[1], a[0]-b[0]))
        return ang + 360 if ang < 0 else ang

    def three_points_intersection(neighbor_vertices, vertices, lines):
        # This function gets a triangle ceter seed which has three neighnors and returns a
        # small mesh triangle around the triangle centers
        mid_point = vertices[neighbor_vertices[0]][1]
        side_point_1 = vertices[neighbor_vertices[1]][1]
        side_point_2 = vertices[neighbor_vertices[2]][1]
        side_point_3 = vertices[neighbor_vertices[3]][1]
        side_points = [side_point_1, side_point_2, side_point_3]

        mid_point_num = neighbor_vertices[0]
        side_point_1_num = neighbor_vertices[1]
        side_point_2_num = neighbor_vertices[2]
        side_point_3_num = neighbor_vertices[3]


        m = np.zeros(3)
        c = np.zeros(3)
        c_up = np.zeros(3)
        c_down = np.zeros(3)
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
            mid_point_12 = lines_Intersect(m[0], c_up[0], m[1], c_down[1])
            mid_point_23 = lines_Intersect(m[1], c_up[1], m[2], c_down[2])
            mid_point_31 = lines_Intersect(m[2], c_up[2], m[0], c_down[0])

        # If counter-clockwise:
        else:
            mid_point_12 = lines_Intersect(m[0], c_down[0], m[1], c_up[1])
            mid_point_23 = lines_Intersect(m[1], c_down[1], m[2], c_up[2])
            mid_point_31 = lines_Intersect(m[2], c_down[2], m[0], c_up[0])


        lines.append([mid_point_num, [side_point_1_num , side_point_2_num, mid_point_12], [side_point_2_num, side_point_3_num, mid_point_23], [side_point_3_num, side_point_1_num, mid_point_31]])
        return lines

#### ###########################################

#  domain = 10


    file_name = "point_seeds"
    open_file = open(file_name, "rb")
    points = pickle.load(open_file)
    open_file.close()

    if (seeds_remove is True): 
        os.remove("point_seeds")  # Deleting the seeds file, created in seeds.py

    vor = Voronoi(points)
    # fig = voronoi_plot_2d(vor)
    # plt.savefig('Voronoi.jpg')
    num = len(points)

    period_neighbor_points_1 = np.zeros((num,2))
    period_neighbor_points_2 = np.zeros((num,2))
    period_neighbor_points_3 = np.zeros((num,2))
    period_neighbor_points_4 = np.zeros((num,2))
    period_neighbor_points_5 = np.zeros((num,2))
    period_neighbor_points_6 = np.zeros((num,2))
    period_neighbor_points_7 = np.zeros((num,2))
    period_neighbor_points_8 = np.zeros((num,2))


    domain_x = domain_y * np.sqrt(3)/1.5/2
    

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

    
    points = np.concatenate((points, period_neighbor_points_1, period_neighbor_points_2, period_neighbor_points_3, period_neighbor_points_4, period_neighbor_points_5, period_neighbor_points_6, period_neighbor_points_7, period_neighbor_points_8))

    tri = Delaunay(points)


    plt.plot(points[:,0], points[:,1], 'o')
    plt.triplot(points[:,0], points[:,1], tri.simplices)


    # Put each triangle corners in one group

    triangles_cor = tri.simplices
    number_of_triangels = len(triangles_cor)

    # Find vertices for each triangle
    # Put each triangle corners and its vertice in one group and make a list, please
    vertices = []
    for i in range(number_of_triangels):
        vertices.append([triangles_cor[i], vertice_generator(points[triangles_cor[i, 0]], points[triangles_cor[i, 1]], points[triangles_cor[i, 2]])])


    # Find the neighbour triangles
    neighbor_triangles_by_number = []
    for i in range(number_of_triangels):
        neighbor_triangles_by_number.append(find_neighbour_triangles(i, triangles_cor[i], triangles_cor))


    # Put neighbour triangles vertices in one group which shows the ridge of Voronoi tessellation

    lines = []      # lines contains each triangle (vertice) number with the side points
    for i in range (len(neighbor_triangles_by_number)):
        if len(neighbor_triangles_by_number[i]) == 2:
            lines = one_points_intersection(neighbor_triangles_by_number[i], vertices, lines)
        if len(neighbor_triangles_by_number[i]) == 3:
            lines = two_points_intersection(neighbor_triangles_by_number[i], vertices, lines)
        if len(neighbor_triangles_by_number[i]) == 4:
            lines = three_points_intersection(neighbor_triangles_by_number[i], vertices, lines)



    surfaces = []

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

                    l = l + 1
                    if point[0] > neighbor[0]:
                        occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(p_up, p_down), occ.addLine(p_down, n_down), occ.addLine(n_down, n_up), occ.addLine(n_up, p_up)])], l)

                    else:
                        occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(p_down, p_up), occ.addLine(p_up, n_up), occ.addLine(n_up, n_down), occ.addLine(n_down, p_down)])], l)

                    surfaces.append([i, neighbor_num, l])

    # shift_y = 0
    # shift_y = -0.2
    # gmsh.initialize()
    occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(occ.addPoint(0, 0+shift_y, 0, lcar), occ.addPoint(domain_x, 0+shift_y, 0, lcar)),\
                                        occ.addLine(occ.addPoint(domain_x, 0+shift_y, 0, lcar), occ.addPoint(domain_x, domain_y+shift_y, 0, lcar)),\
                                        occ.addLine(occ.addPoint(domain_x, domain_y+shift_y, 0, lcar), occ.addPoint(0, domain_y+shift_y, 0, lcar)),\
                                        occ.addLine(occ.addPoint(0, domain_y+shift_y, 0, lcar), occ.addPoint(0, 0+shift_y, 0, lcar))])], 20000)

    base = 'occ.fuse([(2,1)], [(2,2)'
    for i in range(l-2):
        base = base + ',(2,' + str(i+3) + ')'

    command = base + '], 30000)'

    exec(command)

    frame = occ.cut([(2, 20000)], [(2, 30000)])

    

    occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(occ.addPoint(0, 0+shift_y, 0, lcar), occ.addPoint(domain_x, 0+shift_y, 0, lcar)),\
                                        occ.addLine(occ.addPoint(domain_x, 0+shift_y, 0, lcar), occ.addPoint(domain_x, domain_y+shift_y, 0, lcar)),\
                                        occ.addLine(occ.addPoint(domain_x, domain_y+shift_y, 0, lcar), occ.addPoint(0, domain_y+shift_y, 0, lcar)),\
                                        occ.addLine(occ.addPoint(0, domain_y+shift_y, 0, lcar), occ.addPoint(0, 0+shift_y, 0, lcar))])], 20001)

    base = 'occ.cut([(2, 20001)], [(2, 1)'
    for i in range(len(frame[0])-1):
        base = base + ',(2,' + str(i+2) + ')'

    command = base + '], 40000)'
    exec(command)
    gmsh.option.setNumber("Mesh.MeshSizeMax", lcar)

    
    occ.synchronize()
    mat = gmsh.model.addPhysicalGroup(2, [40000])
    # gmsh.model.addPhysicalGroup(2, [40000])
    # gmsh.model.addPhysicalGroup(1, [124], 1)
    # gmsh.model.addPhysicalGroup(1, [169, 176], 2)
    # gmsh.model.addPhysicalGroup(1, [122], 3)
    # gmsh.model.addPhysicalGroup(1, [123], 4)
    # gmsh.model.addPhysicalGroup(1, [155], 5)
    # gmsh.model.addPhysicalGroup(1, [168], 6)
    # gmsh.model.addPhysicalGroup(1, [161, 162], 7)
    # gmsh.model.addPhysicalGroup(1, [164, 165, 166], 8)
    # gmsh.model.addPhysicalGroup(1, [157, 158, 159], 9)
    # gmsh.model.addPhysicalGroup(1, [170, 171, 172, 173, 174, 175], 10)
    # gmsh.model.geo.rotate(gmsh.model.occ.getEntities(40000), 0, 0, 0, 0, 1, 1, math.pi / 4)
    
    zmin = 0; zmax = 0
    ymin = 0+shift_y; ymax = domain_y+shift_y
    xmin = 0; xmax = domain_x
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
                    gmsh.model.mesh.setPeriodic(1, [j[1]], [i[1]], [1, 0, 0, dx,\
                                                                    0, 1, 0, dy,\
                                                                    0, 0, 1, dz,\
                                                                    0, 0, 0, 1 ])

    setPeriodic(0)
    setPeriodic(1)
    # gmsh.fltk.run()

    gmsh.model.mesh.generate()

    gmsh.write(fname + '.msh')
    gmsh.write(fname+"-mesh.vtk")
    gmsh.finalize()

    os.system("gmsh -2 -o " + fname + ".msh -format msh22 " + fname + ".msh")
    os.system("dolfin-convert " + fname + ".msh " +  fname + ".xml")

    mesh = meshio.read(fname+"-mesh.vtk")
    mesh.points = mesh.points[:, :2]
    meshio.write(fname+"-mesh.xdmf", mesh)

    # gmsh.write(fname+"-mesh.vtk")
    # gmsh.finalize()

    # mesh = meshio.read(fname+"-mesh.vtk")
    # mesh.points = mesh.points[:, :2]
    # meshio.write(fname+"-mesh.xdmf", mesh)

    colomn = 2 * row
    number_of_cells = 2 * row * colomn

    mesh = dolfin.Mesh()
    dolfin.XDMFFile(fname+"-mesh.xdmf").read(mesh)
    coord = mesh.coordinates()
    xmax = max(coord[:,0]); xmin = min(coord[:,0])
    ymax = max(coord[:,1]); ymin = min(coord[:,1])
    vol = (xmax - xmin)*(ymax - ymin)
    dV = dolfin.Measure("dx",domain=mesh)
    Vs0 = dolfin.assemble(dolfin.Constant(1) * dV)

    phi = (vol - Vs0)/vol
    print("porosity:" +str(phi))
    return(phi)



def voronoi_tessellation(fname, h, row, domain_y, shift_y, seeds_remove = True):

    import math
    import meshio
    import gmsh
    import os
    import pickle
    import dolfin
    import numpy as np
    from scipy.spatial import Delaunay
    import matplotlib.pyplot as plt
    from scipy.spatial import Voronoi, voronoi_plot_2d


    ########## Initialization#######################################################


    gmsh.initialize()
    gmsh.clear()
    model = gmsh.model
    occ = model.occ

    ############# Parameters #######################################################

    lcar = h/5
    l = 0
    ############# Funtions #########################################################

    def line_btw_points(P, Q):
    #Retunrs the line parameter between two points
        # by = ax + c
        a = Q[1] - P[1]
        b = Q[0] - P[0]
        c = b * P[1] - a * P[0]
        # beta = math.atan(a/b)
        return a, b, c

    def line_2_points(P, Q):
    #Retunrs the line parameter between two points
        m = (Q[1] - P[1])/(Q[0] - P[0])
        c = P[1] - m * P[0]
        return m, c

    def perpendicular_line(P, Q, a, b, c):
    #Returns the perpendicular line to a line passing two points
        mid_point = [(P[0] + Q[0])/2, (P[1] + Q[1])/2]

        # ay = -bx + d
        c = b * (mid_point[0]) + a * (mid_point[1])
        temp = a
        a = -b
        b = temp
        return a, b, c

    def lines_Intersect(m1, c1, m2, c2):
    # Returns the intersection point of two lines
        x = (c2 - c1)/(m1 - m2)
        y = (m2*c1 - m1*c2)/(m2 - m1)
        return [x, y]

    def lines_Intersection(a1, b1, c1, a2, b2, c2):
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


    def vertice_generator(P, Q, S):
        a1, b1, c1 = line_btw_points(P, Q)
        a2, b2, c2 = line_btw_points(Q, S)
        a1, b1, c1 = perpendicular_line(P, Q, a1,b1,c1)
        a2, b2, c2 = perpendicular_line(Q, S, a2,b2,c2)
        x, y = lines_Intersection(a1,b1,c1,a2,b2,c2)
        return [x, y]

    def find_neighbour_triangles(triangle_num, triangle, triangles_cor):
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




        mid_point_up = lines_Intersect(m1, c1_up, m2, c2_down)
        mid_point_down = lines_Intersect(m1, c1_down, m2, c2_up)

        lines.append([mid_point_num, mid_point_up, mid_point_down])
        return lines

    def getAngle(a, b, c):
        ang = math.degrees(math.atan2(c[1]-b[1], c[0]-b[0]) - math.atan2(a[1]-b[1], a[0]-b[0]))
        return ang + 360 if ang < 0 else ang

    def three_points_intersection(neighbor_vertices, vertices, lines):
        # This function gets a triangle ceter seed which has three neighnors and returns a
        # small mesh triangle around the triangle centers
        mid_point = vertices[neighbor_vertices[0]][1]
        side_point_1 = vertices[neighbor_vertices[1]][1]
        side_point_2 = vertices[neighbor_vertices[2]][1]
        side_point_3 = vertices[neighbor_vertices[3]][1]
        side_points = [side_point_1, side_point_2, side_point_3]

        mid_point_num = neighbor_vertices[0]
        side_point_1_num = neighbor_vertices[1]
        side_point_2_num = neighbor_vertices[2]
        side_point_3_num = neighbor_vertices[3]


        m = np.zeros(3)
        c = np.zeros(3)
        c_up = np.zeros(3)
        c_down = np.zeros(3)
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
            mid_point_12 = lines_Intersect(m[0], c_up[0], m[1], c_down[1])
            mid_point_23 = lines_Intersect(m[1], c_up[1], m[2], c_down[2])
            mid_point_31 = lines_Intersect(m[2], c_up[2], m[0], c_down[0])

        # If counter-clockwise:
        else:
            mid_point_12 = lines_Intersect(m[0], c_down[0], m[1], c_up[1])
            mid_point_23 = lines_Intersect(m[1], c_down[1], m[2], c_up[2])
            mid_point_31 = lines_Intersect(m[2], c_down[2], m[0], c_up[0])


        lines.append([mid_point_num, [side_point_1_num , side_point_2_num, mid_point_12], [side_point_2_num, side_point_3_num, mid_point_23], [side_point_3_num, side_point_1_num, mid_point_31]])
        return lines

#### ###########################################

#  domain = 10


    file_name = "point_seeds_random"
    # file_name = "point_seeds"
    open_file = open(file_name, "rb")
    points = pickle.load(open_file)
    open_file.close()

    if (seeds_remove is True): 
        os.remove("point_seeds")  # Deleting the seeds file, created in seeds.py

    vor = Voronoi(points)
    # fig = voronoi_plot_2d(vor)
    # plt.savefig('Voronoi.jpg')
    num = len(points)

    period_neighbor_points_1 = np.zeros((num,2))
    period_neighbor_points_2 = np.zeros((num,2))
    period_neighbor_points_3 = np.zeros((num,2))
    period_neighbor_points_4 = np.zeros((num,2))
    period_neighbor_points_5 = np.zeros((num,2))
    period_neighbor_points_6 = np.zeros((num,2))
    period_neighbor_points_7 = np.zeros((num,2))
    period_neighbor_points_8 = np.zeros((num,2))


    domain_x = domain_y * np.sqrt(3)/1.5
    

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

    
    points = np.concatenate((points, period_neighbor_points_1, period_neighbor_points_2, period_neighbor_points_3, period_neighbor_points_4, period_neighbor_points_5, period_neighbor_points_6, period_neighbor_points_7, period_neighbor_points_8))

    tri = Delaunay(points)


    plt.plot(points[:,0], points[:,1], 'o')
    plt.triplot(points[:,0], points[:,1], tri.simplices)


    # Put each triangle corners in one group

    triangles_cor = tri.simplices
    number_of_triangels = len(triangles_cor)

    # Find vertices for each triangle
    # Put each triangle corners and its vertice in one group and make a list, please
    vertices = []
    for i in range(number_of_triangels):
        vertices.append([triangles_cor[i], vertice_generator(points[triangles_cor[i, 0]], points[triangles_cor[i, 1]], points[triangles_cor[i, 2]])])


    # Find the neighbour triangles
    neighbor_triangles_by_number = []
    for i in range(number_of_triangels):
        neighbor_triangles_by_number.append(find_neighbour_triangles(i, triangles_cor[i], triangles_cor))


    # Put neighbour triangles vertices in one group which shows the ridge of Voronoi tessellation

    lines = []      # lines contains each triangle (vertice) number with the side points
    for i in range (len(neighbor_triangles_by_number)):
        if len(neighbor_triangles_by_number[i]) == 2:
            lines = one_points_intersection(neighbor_triangles_by_number[i], vertices, lines)
        if len(neighbor_triangles_by_number[i]) == 3:
            lines = two_points_intersection(neighbor_triangles_by_number[i], vertices, lines)
        if len(neighbor_triangles_by_number[i]) == 4:
            lines = three_points_intersection(neighbor_triangles_by_number[i], vertices, lines)



    surfaces = []

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

                    l = l + 1
                    if point[0] > neighbor[0]:
                        occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(p_up, p_down), occ.addLine(p_down, n_down), occ.addLine(n_down, n_up), occ.addLine(n_up, p_up)])], l)

                    else:
                        occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(p_down, p_up), occ.addLine(p_up, n_up), occ.addLine(n_up, n_down), occ.addLine(n_down, p_down)])], l)

                    surfaces.append([i, neighbor_num, l])

    # shift_y = 0
    # shift_y = -0.2
    # gmsh.initialize()
    occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(occ.addPoint(0, 0+shift_y, 0, lcar), occ.addPoint(domain_x, 0+shift_y, 0, lcar)),\
                                        occ.addLine(occ.addPoint(domain_x, 0+shift_y, 0, lcar), occ.addPoint(domain_x, domain_y+shift_y, 0, lcar)),\
                                        occ.addLine(occ.addPoint(domain_x, domain_y+shift_y, 0, lcar), occ.addPoint(0, domain_y+shift_y, 0, lcar)),\
                                        occ.addLine(occ.addPoint(0, domain_y+shift_y, 0, lcar), occ.addPoint(0, 0+shift_y, 0, lcar))])], 20000)

    base = 'occ.fuse([(2,1)], [(2,2)'
    for i in range(l-2):
        base = base + ',(2,' + str(i+3) + ')'

    # base = 'occ.cut([(2,20000)], [(2,2)'
    # for i in range(l-2):
    #     base = base + ',(2,' + str(i+3) + ')'
    
    command = base + '], 30000)'

    exec(command)

    frame = occ.cut([(2, 20000)], [(2, 30000)])

    # occ.synchronize()
    # gmsh.fltk.run()
    

    occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(occ.addPoint(0, 0+shift_y, 0, lcar), occ.addPoint(domain_x, 0+shift_y, 0, lcar)),\
                                        occ.addLine(occ.addPoint(domain_x, 0+shift_y, 0, lcar), occ.addPoint(domain_x, domain_y+shift_y, 0, lcar)),\
                                        occ.addLine(occ.addPoint(domain_x, domain_y+shift_y, 0, lcar), occ.addPoint(0, domain_y+shift_y, 0, lcar)),\
                                        occ.addLine(occ.addPoint(0, domain_y+shift_y, 0, lcar), occ.addPoint(0, 0+shift_y, 0, lcar))])], 20001)

    base = 'occ.cut([(2, 20001)], [(2, 1)'
    for i in range(len(frame[0])-1):
        base = base + ',(2,' + str(i+2) + ')'

    command = base + '], 40000)'
    exec(command)
    
    gmsh.option.setNumber("Mesh.MeshSizeMax", lcar)

    
    occ.synchronize()
    mat = gmsh.model.addPhysicalGroup(2, [40000])
    # gmsh.model.geo.rotate(gmsh.model.occ.getEntities(40000), 0, 0, 0, 0, 1, 1, math.pi / 4)
    
    zmin = 0; zmax = 0
    ymin = 0+shift_y; ymax = domain_y+shift_y
    xmin = 0; xmax = domain_x
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
                    gmsh.model.mesh.setPeriodic(1, [j[1]], [i[1]], [1, 0, 0, dx,\
                                                                    0, 1, 0, dy,\
                                                                    0, 0, 1, dz,\
                                                                    0, 0, 0, 1 ])

    setPeriodic(0)
    setPeriodic(1)
    gmsh.model.mesh.generate()

    gmsh.write(fname + '.msh')
    gmsh.write(fname+"-mesh.vtk")
    gmsh.finalize()

    os.system("gmsh -2 -o " + fname + ".msh -format msh22 " + fname + ".msh")
    os.system("dolfin-convert " + fname + ".msh " +  fname + ".xml")

    mesh = meshio.read(fname+"-mesh.vtk")
    mesh.points = mesh.points[:, :2]
    meshio.write(fname+"-mesh.xdmf", mesh)

    # gmsh.write(fname+"-mesh.vtk")
    # gmsh.finalize()

    # mesh = meshio.read(fname+"-mesh.vtk")
    # mesh.points = mesh.points[:, :2]
    # meshio.write(fname+"-mesh.xdmf", mesh)

    colomn = 2 * row
    number_of_cells = 2 * row * colomn

    mesh = dolfin.Mesh()
    dolfin.XDMFFile(fname+"-mesh.xdmf").read(mesh)
    coord = mesh.coordinates()
    xmax = max(coord[:,0]); xmin = min(coord[:,0])
    ymax = max(coord[:,1]); ymin = min(coord[:,1])
    vol = (xmax - xmin)*(ymax - ymin)
    dV = dolfin.Measure("dx",domain=mesh)
    Vs0 = dolfin.assemble(dolfin.Constant(1) * dV)

    phi = (vol - Vs0)/vol
    print("porosity:" +str(phi))
    return(phi)


def Hallow_Square(fname, width, r, shift_x, shift_y):

    import gmsh
    import os
    import math
    import meshio
    import dolfin

    ###################################################################################################################################################################
    center = width/2
    xmin = 0.
    ymin = 0.
    zmin = 0.
    xmax = width
    ymax = width
    zmax = width
    x0 = center
    y0 = center
    z0 = center
    r0 = r
    l = width/50 
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
    hole_tag = 2
    hole_tag2 = 3
    hole_tag3 = 4
    hole_tag4 = 5
    hole_tag5 = 6
    rve_tag = 7

    gmsh.model.occ.addRectangle(x=xmin+shift_x, y=ymin+shift_y, z=0, dx=xmax-xmin, dy=ymax-ymin, tag=box_tag)
    # gmsh.model.occ.addDisk(xc=x0, yc=y0, zc=0, rx=r0, ry=r0, tag=hole_tag)
    gmsh.model.occ.addDisk(xc=xmin, yc=ymin, zc=0, rx=r0, ry=r0, tag=hole_tag2)
    gmsh.model.occ.addDisk(xc=xmax, yc=ymin, zc=0, rx=r0, ry=r0, tag=hole_tag3)
    gmsh.model.occ.addDisk(xc=xmax, yc=ymax, zc=0, rx=r0, ry=r0, tag=hole_tag4)
    gmsh.model.occ.addDisk(xc=xmin, yc=ymax, zc=0, rx=r0, ry=r0, tag=hole_tag5)
    # gmsh.model.occ.cut(objectDimTags=[(2, box_tag)], toolDimTags=[(2, hole_tag)], tag=rve_tag)
    # gmsh.model.occ.cut(objectDimTags=[(2, box_tag)], toolDimTags=[(2, hole_tag), (2, hole_tag2)], tag=rve_tag)
    # gmsh.model.occ.cut(objectDimTags=[(2, box_tag)], toolDimTags=[(2, hole_tag), (2, hole_tag2), (2, hole_tag3), (2, hole_tag4), (2, hole_tag5)], tag=rve_tag)
    gmsh.model.occ.cut(objectDimTags=[(2, box_tag)], toolDimTags=[(2, hole_tag2), (2, hole_tag3), (2, hole_tag4), (2, hole_tag5)], tag=rve_tag)
    # gmsh.model.geo.rotate(rve_tag, 0, 0, 0, 0, 0, 1, math.pi / 4)
    gmsh.model.occ.synchronize()
    gmsh.model.addPhysicalGroup(dim=2, tags=[rve_tag])
    # gmsh.model.mesh.setPeriodic(dim=1, tags=[2], tagsMaster=[1], affineTransform=[1, 0, 0, xmax-xmin,\
    #                                                                                 0, 1, 0, 0        ,\
    #                                                                                 0, 0, 1, 0        ,\
    #                                                                                 0, 0, 0, 1        ])
    # gmsh.model.mesh.setPeriodic(dim=1, tags=[4], tagsMaster=[3], affineTransform=[1, 0, 0, 0        ,\
    #                                                                                 0, 1, 0, ymax-ymin,\
    #                                                                                 0, 0, 1, 0        ,\
    #                                                                                 0, 0, 0, 1        ])

    
    setPeriodic(0)
    setPeriodic(1)
    gmsh.model.mesh.setSize(dimTags=gmsh.model.getEntities(0), size=l)
    gmsh.model.mesh.generate(dim=2)


    # gmsh.model.mesh.generate()
    gmsh.write(fname + '.msh')
    os.system("gmsh -2 -o " + fname + ".msh -format msh22 " + fname + ".msh")
    os.system("dolfin-convert " + fname + ".msh " +  fname + ".xml")

    gmsh.write(fname+"-mesh.vtk")
    gmsh.finalize()

    mesh = meshio.read(fname+"-mesh.vtk")
    mesh.points = mesh.points[:, :2]
    meshio.write(fname+"-mesh.xdmf", mesh)

    mesh = dolfin.Mesh()
    dolfin.XDMFFile(fname+"-mesh.xdmf").read(mesh)
    coord = mesh.coordinates()
    xmax = max(coord[:,0]); xmin = min(coord[:,0])
    ymax = max(coord[:,1]); ymin = min(coord[:,1])
    vol = (xmax - xmin)*(ymax - ymin)
    dV = dolfin.Measure("dx",domain=mesh)
    Vs0 = dolfin.assemble(dolfin.Constant(1) * dV)

    phi = (vol - Vs0)/vol
    print("porosity:" +str(phi))

def Porous_media(fname, width, r, shift_x, shift_y):

    import gmsh
    import os
    import math
    import meshio
    import dolfin

    ###################################################################################################################################################################
    center = width/2
    xmin = 0.
    ymin = 0.
    zmin = 0.
    xmax = width
    ymax = width
    zmax = width
    x0 = center
    y0 = center
    z0 = center
    r0 = r
    l = width/50 
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
    hole_tag = 2
    hole_tag2 = 3
    hole_tag3 = 4
    hole_tag4 = 5
    hole_tag5 = 6
    rve_tag = 7

    gmsh.model.occ.addRectangle(x=xmin+shift_x, y=ymin+shift_y, z=0, dx=xmax-xmin, dy=ymax-ymin, tag=box_tag)
    # gmsh.model.occ.addDisk(xc=x0, yc=y0, zc=0, rx=r0, ry=r0, tag=hole_tag)
    gmsh.model.occ.addDisk(xc=xmin, yc=(ymin+ymax)*5/8, zc=0, rx=r0, ry=r0, tag=hole_tag2)
    gmsh.model.occ.addDisk(xc=xmin, yc=ymin, zc=0, rx=r0, ry=r0, tag=hole_tag3)
    gmsh.model.occ.addDisk(xc=(xmax + xmin)*5/8, yc=(ymax + ymin)/2, zc=0, rx=r0, ry=r0, tag=hole_tag4)
    gmsh.model.occ.addDisk(xc=xmax, yc=(ymax+ymin)*3/4, zc=0, rx=r0, ry=r0, tag=hole_tag5)
    # gmsh.model.occ.cut(objectDimTags=[(2, box_tag)], toolDimTags=[(2, hole_tag)], tag=rve_tag)
    # gmsh.model.occ.cut(objectDimTags=[(2, box_tag)], toolDimTags=[(2, hole_tag), (2, hole_tag2)], tag=rve_tag)
    # gmsh.model.occ.cut(objectDimTags=[(2, box_tag)], toolDimTags=[(2, hole_tag), (2, hole_tag2), (2, hole_tag3), (2, hole_tag4), (2, hole_tag5)], tag=rve_tag)
    gmsh.model.occ.cut(objectDimTags=[(2, box_tag)], toolDimTags=[(2, hole_tag2), (2, hole_tag3), (2, hole_tag4), (2, hole_tag5)], tag=rve_tag)
    # gmsh.model.geo.rotate(rve_tag, 0, 0, 0, 0, 0, 1, math.pi / 4)
    gmsh.model.occ.synchronize()
    gmsh.model.addPhysicalGroup(dim=2, tags=[rve_tag])
    # gmsh.model.mesh.setPeriodic(dim=1, tags=[2], tagsMaster=[1], affineTransform=[1, 0, 0, xmax-xmin,\
    #                                                                                 0, 1, 0, 0        ,\
    #                                                                                 0, 0, 1, 0        ,\
    #                                                                                 0, 0, 0, 1        ])
    # gmsh.model.mesh.setPeriodic(dim=1, tags=[4], tagsMaster=[3], affineTransform=[1, 0, 0, 0        ,\
    #                                                                                 0, 1, 0, ymax-ymin,\
    #                                                                                 0, 0, 1, 0        ,\
    #                                                                                 0, 0, 0, 1        ])

    
    setPeriodic(0)
    setPeriodic(1)
    gmsh.model.mesh.setSize(dimTags=gmsh.model.getEntities(0), size=l)
    gmsh.model.mesh.generate(dim=2)


    # gmsh.model.mesh.generate()
    gmsh.write(fname + '.msh')
    os.system("gmsh -2 -o " + fname + ".msh -format msh22 " + fname + ".msh")
    os.system("dolfin-convert " + fname + ".msh " +  fname + ".xml")

    gmsh.write(fname+"-mesh.vtk")
    gmsh.finalize()

    mesh = meshio.read(fname+"-mesh.vtk")
    mesh.points = mesh.points[:, :2]
    meshio.write(fname+"-mesh.xdmf", mesh)

    mesh = dolfin.Mesh()
    dolfin.XDMFFile(fname+"-mesh.xdmf").read(mesh)
    coord = mesh.coordinates()
    xmax = max(coord[:,0]); xmin = min(coord[:,0])
    ymax = max(coord[:,1]); ymin = min(coord[:,1])
    vol = (xmax - xmin)*(ymax - ymin)
    dV = dolfin.Measure("dx",domain=mesh)
    Vs0 = dolfin.assemble(dolfin.Constant(1) * dV)

    phi = (vol - Vs0)/vol
    print("porosity:" +str(phi))

def Hallow_cube(fname, width, r, shift_x, shift_y, shift_z):

    import gmsh
    import os
    import math
    import meshio
    import dolfin

    ###################################################################################################################################################################
    center = width/2
    xmin = 0.
    ymin = 0.
    zmin = 0.
    xmax = width
    ymax = width
    zmax = width
    x0 = center
    y0 = center
    z0 = center
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
    gmsh.write(fname + '.msh')
    os.system("gmsh -2 -o " + fname + ".msh -format msh22 " + fname + ".msh")
    os.system("dolfin-convert " + fname + ".msh " +  fname + ".xml")

    gmsh.write(fname+"-mesh.vtk")
    gmsh.finalize()

    mesh = meshio.read(fname+"-mesh.vtk")
    mesh.points = mesh.points[:, :3]
    meshio.write(fname+"-mesh.xdmf", mesh)

    mesh = dolfin.Mesh()
    dolfin.XDMFFile(fname+"-mesh.xdmf").read(mesh)
    
    vol = width**3
    phi = (4/3*math.pi*r**3)/vol
    print("porosity:" +str(phi))

def Hallow_array(fname, shift_x, shift_y, width):

    import gmsh
    import os
    import math
    import meshio
    import dolfin

    ###################################################################################################################################################################

    # dist = 0.25
    # width = 3 * dist
    dist = width/3
    # shift_x = dist/2
    # shift_y = dist/2

    xmin = 0.
    ymin = 0.
    zmin = 0.
    xmax = width
    ymax = width
    zmax = width
    r0 = dist/5
    l = dist/50
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
    hole_tag11 = 2
    hole_tag12 = 3
    hole_tag13 = 4
    hole_tag14 = 5

    hole_tag21 = 6
    hole_tag22 = 7
    hole_tag23 = 8
    hole_tag24 = 9

    hole_tag31 = 10
    hole_tag32 = 11
    hole_tag33 = 12
    hole_tag34 = 13

    hole_tag41 = 14
    hole_tag42 = 15
    hole_tag43 = 16
    hole_tag44 = 17

    rve_tag = 18

    gmsh.model.occ.addRectangle(x=xmin + shift_x, y=ymin + shift_y, z=0, dx=xmax-xmin, dy=ymax-ymin, tag=box_tag)

    gmsh.model.occ.addDisk(xc=xmin, yc=ymin, zc=0, rx=r0, ry=r0, tag=hole_tag11)
    gmsh.model.occ.addDisk(xc=xmin + dist, yc=ymin, zc=0, rx=r0, ry=r0, tag=hole_tag12)
    gmsh.model.occ.addDisk(xc=xmin + 2*dist, yc=ymin, zc=0, rx=r0, ry=r0, tag=hole_tag13)
    gmsh.model.occ.addDisk(xc=xmax, yc=ymin, zc=0, rx=r0, ry=r0, tag=hole_tag14)

    gmsh.model.occ.addDisk(xc=xmin, yc=ymin + dist, zc=0, rx=r0, ry=r0, tag=hole_tag21)
    gmsh.model.occ.addDisk(xc=xmin + dist, yc=ymin + dist, zc=0, rx=r0, ry=r0, tag=hole_tag22)
    gmsh.model.occ.addDisk(xc=xmin + 2*dist, yc=ymin + dist, zc=0, rx=r0, ry=r0, tag=hole_tag23)
    gmsh.model.occ.addDisk(xc=xmax, yc=ymin + dist, zc=0, rx=r0, ry=r0, tag=hole_tag24)

    gmsh.model.occ.addDisk(xc=xmin, yc=ymin + 2* dist, zc=0, rx=r0, ry=r0, tag=hole_tag31)
    gmsh.model.occ.addDisk(xc=xmin + dist, yc=ymin + 2*dist, zc=0, rx=r0, ry=r0, tag=hole_tag32)
    gmsh.model.occ.addDisk(xc=xmin + 2*dist, yc=ymin + 2*dist, zc=0, rx=r0, ry=r0, tag=hole_tag33)
    gmsh.model.occ.addDisk(xc=xmax, yc=ymin + 2*dist, zc=0, rx=r0, ry=r0, tag=hole_tag34)

    gmsh.model.occ.addDisk(xc=xmin, yc=ymax, zc=0, rx=r0, ry=r0, tag=hole_tag41)
    gmsh.model.occ.addDisk(xc=xmin + dist, yc=ymax, zc=0, rx=r0, ry=r0, tag=hole_tag42)
    gmsh.model.occ.addDisk(xc=xmin + 2*dist, yc=ymax, zc=0, rx=r0, ry=r0, tag=hole_tag43)
    gmsh.model.occ.addDisk(xc=xmax, yc=ymax, zc=0, rx=r0, ry=r0, tag=hole_tag44)

    

    
    # gmsh.model.occ.cut(objectDimTags=[(2, box_tag)], toolDimTags=[(2, hole_tag)], tag=rve_tag)
    # gmsh.model.occ.cut(objectDimTags=[(2, box_tag)], toolDimTags=[(2, hole_tag), (2, hole_tag2)], tag=rve_tag)
    # gmsh.model.occ.cut(objectDimTags=[(2, box_tag)], toolDimTags=[(2, hole_tag), (2, hole_tag2), (2, hole_tag3), (2, hole_tag4), (2, hole_tag5)], tag=rve_tag)
    gmsh.model.occ.cut(objectDimTags=[(2, box_tag)], toolDimTags=[(2, hole_tag11), (2, hole_tag12), (2, hole_tag13), (2, hole_tag14),\
                                                                  (2, hole_tag21), (2, hole_tag22), (2, hole_tag23), (2, hole_tag24),
                                                                  (2, hole_tag31), (2, hole_tag32), (2, hole_tag33), (2, hole_tag34),
                                                                  (2, hole_tag41), (2, hole_tag42), (2, hole_tag43), (2, hole_tag44)], tag=rve_tag)
    # gmsh.model.geo.rotate(rve_tag, 0, 0, 0, 0, 0, 1, math.pi / 4)
    gmsh.model.occ.synchronize()
    gmsh.model.addPhysicalGroup(dim=2, tags=[rve_tag])
    # gmsh.model.mesh.setPeriodic(dim=1, tags=[2], tagsMaster=[1], affineTransform=[1, 0, 0, xmax-xmin,\
    #                                                                                 0, 1, 0, 0        ,\
    #                                                                                 0, 0, 1, 0        ,\
    #                                                                                 0, 0, 0, 1        ])
    # gmsh.model.mesh.setPeriodic(dim=1, tags=[4], tagsMaster=[3], affineTransform=[1, 0, 0, 0        ,\
    #                                                                                 0, 1, 0, ymax-ymin,\
    #                                                                                 0, 0, 1, 0        ,\
    #                                                                                 0, 0, 0, 1        ])

    
    setPeriodic(0)
    setPeriodic(1)
    gmsh.model.mesh.setSize(dimTags=gmsh.model.getEntities(0), size=l)
    gmsh.model.mesh.generate(dim=2)


    # gmsh.model.mesh.generate()
    gmsh.write(fname + '.msh')
    os.system("gmsh -2 -o " + fname + ".msh -format msh22 " + fname + ".msh")
    os.system("dolfin-convert " + fname + ".msh " +  fname + ".xml")

    gmsh.write(fname+"-mesh.vtk")
    gmsh.finalize()

    mesh = meshio.read(fname+"-mesh.vtk")
    mesh.points = mesh.points[:, :2]
    meshio.write(fname+"-mesh.xdmf", mesh)




def Hexag_incl(fname, width, angle, r, shift_x, shift_y):

    import math
    import gmsh
    import os
    import meshio

    gmsh.initialize()
    gmsh.clear()
    model = gmsh.model
    occ = model.occ

        ############# Parameters #######################################################
    a = width
    R = r
    t = angle
    b = math.sin(t)*a
    c = a * math.cos(t)

    lcar = a/50
        ############# Funtions #########################################################
    xmin = 0; xmax = c
    ymin = 0; ymax = b
    zmin = 0; zmax = 0
    e = 1e-6

    occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(occ.addPoint(0, 0, 0, lcar), occ.addPoint(a, 0, 0, lcar)),\
                                        occ.addLine(occ.addPoint(a, 0, 0, lcar), occ.addPoint(a+c, b, 0, lcar)),\
                                        occ.addLine(occ.addPoint(a+c, b, 0, lcar), occ.addPoint(c, b, 0, lcar)),\
                                        occ.addLine(occ.addPoint(c, b, 0, lcar), occ.addPoint(0, 0, 0, lcar))])], 1)

    # occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(occ.addPoint(0+3*r, 0+3*r, 0, lcar), occ.addPoint(a+3*r, 0+3*r, 0, lcar)),\
    #                                     occ.addLine(occ.addPoint(a+3*r, 0+3*r, 0, lcar), occ.addPoint(a+c+3*r, b+3*r, 0, lcar)),\
    #                                     occ.addLine(occ.addPoint(a+c+3*r, b+3*r, 0, lcar), occ.addPoint(c+3*r, b+3*r, 0, lcar)),\
    #                                     occ.addLine(occ.addPoint(c+3*r, b+3*r, 0, lcar), occ.addPoint(0+3*r, 0+3*r, 0, lcar))])], 1)
    
    occ.addDisk(xc=0+shift_x, yc=0+shift_y, zc=0, rx=R, ry=R, tag=2)
    occ.addDisk(xc=a+shift_x, yc=0+shift_y, zc=0, rx=R, ry=R, tag=3)
    occ.addDisk(xc=a+c+shift_x, yc=b+shift_y, zc=0, rx=R, ry=R, tag=4)
    occ.addDisk(xc=c+shift_x, yc=b+shift_y, zc=0, rx=R, ry=R, tag=5)
    occ.cut(objectDimTags=[(2, 1)], toolDimTags=[(2, 2), (2, 3), (2, 4), (2, 5)], tag=6)
    occ.synchronize()
    gmsh.model.addPhysicalGroup(2, [6])

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


    setPeriodic(0)
    setPeriodic(1)
    gmsh.model.mesh.setSize(dimTags=gmsh.model.getEntities(0), size=lcar)
    gmsh.model.mesh.generate(dim=2)



    # gmsh.model.mesh.generate()
    gmsh.write(fname + '.msh')
    os.system("gmsh -2 -o " + fname + ".msh -format msh22 " + fname + ".msh")
    os.system("dolfin-convert " + fname + ".msh " +  fname + ".xml")

    gmsh.write(fname+"-mesh.vtk")
    gmsh.finalize()

    mesh = meshio.read(fname+"-mesh.vtk")
    mesh.points = mesh.points[:, :2]
    meshio.write(fname+"-mesh.xdmf", mesh)



def Hexg(fname, phi, R):

    import math
    import gmsh
    import os
    import meshio

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
    gmsh.write(fname + '.msh')
    os.system("gmsh -2 -o " + fname + ".msh -format msh22 " + fname + ".msh")
    os.system("dolfin-convert " + fname + ".msh " +  fname + ".xml")

    gmsh.write(fname+"-mesh.vtk")
    gmsh.finalize()

    mesh = meshio.read(fname+"-mesh.vtk")
    mesh.points = mesh.points[:, :2]
    meshio.write(fname+"-mesh.xdmf", mesh)


def Hex_rectangle(fname, width, r, shift_x, shift_y):
    

    import gmsh
    import os
    import math
    import meshio
    import dolfin

    ###################################################################################################################################################################
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
    l = width/40
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
    # gmsh.model.occ.cut(objectDimTags=[(2, box_tag)], toolDimTags=[(2, hole_tag)], tag=rve_tag)
    # gmsh.model.occ.cut(objectDimTags=[(2, box_tag)], toolDimTags=[(2, hole_tag), (2, hole_tag2)], tag=rve_tag)
    gmsh.model.occ.cut(objectDimTags=[(2, box_tag)], toolDimTags=[(2, hole_tag), (2, hole_tag2), (2, hole_tag3), (2, hole_tag4), (2, hole_tag5), (2, hole_tag6)], tag=rve_tag)
    # gmsh.model.occ.cut(objectDimTags=[(2, box_tag)], toolDimTags=[(2, hole_tag2), (2, hole_tag3), (2, hole_tag4), (2, hole_tag5)], tag=rve_tag)
    # gmsh.model.geo.rotate(rve_tag, 0, 0, 0, 0, 0, 1, math.pi / 4)
    gmsh.model.occ.synchronize()
    gmsh.model.addPhysicalGroup(dim=2, tags=[rve_tag])
    # gmsh.model.mesh.setPeriodic(dim=1, tags=[2], tagsMaster=[1], affineTransform=[1, 0, 0, xmax-xmin,\
    #                                                                                 0, 1, 0, 0        ,\
    #                                                                                 0, 0, 1, 0        ,\
    #                                                                                 0, 0, 0, 1        ])
    # gmsh.model.mesh.setPeriodic(dim=1, tags=[4], tagsMaster=[3], affineTransform=[1, 0, 0, 0        ,\
    #                                                                                 0, 1, 0, ymax-ymin,\
    #                                                                                 0, 0, 1, 0        ,\
    #                                                                                 0, 0, 0, 1        ])

    
    setPeriodic(0)
    setPeriodic(1)
    gmsh.model.mesh.setSize(dimTags=gmsh.model.getEntities(0), size=l)
    # gmsh.model.mesh.generate(dim=2)


    gmsh.model.mesh.generate()
    gmsh.write(fname + '.msh')
    os.system("gmsh -2 -o " + fname + ".msh -format msh22 " + fname + ".msh")
    os.system("dolfin-convert " + fname + ".msh " +  fname + ".xml")

    gmsh.write(fname+"-mesh.vtk")
    gmsh.finalize()

    mesh = meshio.read(fname+"-mesh.vtk")
    mesh.points = mesh.points[:, :2]
    meshio.write(fname+"-mesh.xdmf", mesh)

    mesh = dolfin.Mesh()
    dolfin.XDMFFile(fname+"-mesh.xdmf").read(mesh)
    coord = mesh.coordinates()
    xmax = max(coord[:,0]); xmin = min(coord[:,0])
    ymax = max(coord[:,1]); ymin = min(coord[:,1])
    vol = (xmax - xmin)*(ymax - ymin)
    dV = dolfin.Measure("dx",domain=mesh)
    Vs0 = dolfin.assemble(dolfin.Constant(1) * dV)

    phi = (vol - Vs0)/vol
    print("porosity:" +str(phi))



def TKD(fname, l, h):
    import gmsh
    import os
    import math
    import meshio
    import dolfin

    
    ###################################################################################################################################################################
    l =1
    h = 0.1
    lcar = h/5

    gmsh.initialize()
    gmsh.clear()
    model = gmsh.model
    occ = model.occ
    occ.mesh.ToleranceInitialDelaunay=1e-12
    model.mesh.ToleranceInitialDelaunay=1e-12

    # occ.mesh.MinimumCircleNodes=16
    # occ.mesh.MinimumCurveNodes=16
    # occ.mesh.MinimumCirclePoints=16
    # occ.mesh.MinimumCurvePoints=16
    # occ.mesh.MeshSizeFromCurvature=16
    # occ.mesh.MeshSizeFromPoints=0
    # occ.mesh.MeshSizeFromParametricPoints=0

    p1 = occ.addPoint(0, 0, l, lcar)
    p2 = occ.addPoint(0, 0, -l, lcar)
    p3 = occ.addPoint(l, 0, 0, lcar)
    p4 = occ.addPoint(0, -l, 0, lcar)
    p5 = occ.addPoint(-l, 0, 0, lcar)
    p6 = occ.addPoint(0, l, 0, lcar)

    s1 = occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(p1, p3), occ.addLine(p3, p6), occ.addLine(p6, p1)])])
    s2 = occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(p1, p6), occ.addLine(p6, p5), occ.addLine(p5, p1)])])
    s3 = occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(p1, p5), occ.addLine(p5, p4), occ.addLine(p4, p1)])])
    s4 = occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(p1, p4), occ.addLine(p4, p3), occ.addLine(p3, p1)])])
    s5 = occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(p3, p6), occ.addLine(p3, p2), occ.addLine(p2, p6)])])
    s6 = occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(p5, p6), occ.addLine(p6, p3), occ.addLine(p3, p5)])])
    s7 = occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(p4, p5), occ.addLine(p5, p2), occ.addLine(p2, p4)])])
    s8 = occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(p4, p3), occ.addLine(p3, p2), occ.addLine(p2, p4)])])

    occ.addVolume([occ.addSurfaceLoop([s1, s2, s3, s4, s5, s6, s7, s8])], 100)

    
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
                    gmsh.model.mesh.setPeriodic(1, [j[1]], [i[1]], [1, 0, 0, dx,\
                                                                    0, 1, 0, dy,\
                                                                    0, 0, 1, dz,\
                                                                    0, 0, 0, 1 ])

    # setPeriodic(0)
    # setPeriodic(1)
    occ.mesh.ToleranceInitialDelaunay=1e-12
    occ.synchronize()
    gmsh.model.mesh.ToleranceInitialDelaunay=1e-12
    gmsh.fltk.run()
    model.mesh.generate()

    gmsh.write(fname + '.msh')
    gmsh.write(fname+"-mesh.vtk")
    gmsh.finalize()

    os.system("gmsh -2 -o " + fname + ".msh -format msh22 " + fname + ".msh")
    os.system("dolfin-convert " + fname + ".msh " +  fname + ".xml")

    mesh = meshio.read(fname+"-mesh.vtk")
    mesh.points = mesh.points[:, :2]
    meshio.write(fname+"-mesh.xdmf", mesh)

