################################################################################
###                                                                          ###
### Created by Mahdi Manoochehrtayebi, 2021-2022                             ###
### Supervisors: Dr. Aline Bel-Brunon, Dr. Martin Genet                      ###
###                                                                          ###
### École Polytechnique, Palaiseau, France                                   ###
###                                                                          ###
################################################################################

import pickle
import numpy as np
################################################################################

def semi_regular(DoI, row, domain_y):    

    import numpy as np
    import random

    def odd_position(colomn_number, row_number):
        #This funcrtion creates points in a row for odd colomns

        cell_center = [cell_width * colomn_number + cell_width/2,  (2 * cell_height * row_number + cell_height/2)]

        left_treshold = cell_center[0] - width_treshold/2
        lower_treshold = cell_center[1] - width_treshold/2

        random_point = [left_treshold + width_treshold * random.random(), lower_treshold + height_treshold * random.random()]
        return(random_point)

    def even_position(colomn_number, row_number):
        #This funcrtion creates points in a row for even colomns
        cell_center = [cell_width * colomn_number, (2* cell_height * row_number + 3*cell_height/2)]

        left_treshold = cell_center[0] - width_treshold/2
        lower_treshold = cell_center[1] - width_treshold/2

        if colomn_number == 1:
            random_point = [cell_center[0] + width_treshold * random.random()/2, lower_treshold + height_treshold * random.random()]

        elif colomn_number == colomn:
            random_point = [left_treshold + width_treshold * random.random()/2, lower_treshold + height_treshold * random.random()]

        else:
            random_point = [left_treshold + width_treshold * random.random(), lower_treshold + height_treshold * random.random()]

        return(random_point)


    domain_x = domain_y * np.sqrt(3)/1.5 
    colomn = 2 * row
    
    points = np.zeros(((2 * colomn ) * row,2))
    cell_height = domain_y/row/2
    cell_width = domain_x/colomn 

    width_treshold = cell_width * DoI
    height_treshold = cell_height * DoI


    k = 0
    for i in range(row):
        for j in range(colomn):

            points[2*k] = odd_position(j, i)
            points[2*k+1] = even_position(j, i)
            k += 1
    print(points)

    # Savig the points to use them for generating the mesh
    file_name = "point_seeds"
    open_file = open(file_name, "wb")
    pickle.dump(points, open_file)
    open_file.close()



# This module generates random seeds
def random(seeds_nom):
    points = np.random.rand(seeds_nom,2)

    file_name = "point_seeds"
    open_file = open(file_name, "wb")
    pickle.dump(points, open_file)
    open_file.close()
