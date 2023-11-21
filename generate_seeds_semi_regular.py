################################################################################
###                                                                          ###
### Created by Mahdi Manoochehrtayebi, 2021-2022                             ###
### Supervisors: Dr. Aline Bel-Brunon, Dr. Martin Genet                      ###
###                                                                          ###
### École Polytechnique, Palaiseau, France                                   ###
###                                                                          ###
################################################################################

import numpy
import random

import micro_poro_structure_generator as gen

################################################################################

def generate_seeds_semi_regular(
        DoI,
        row,
        domain_y,
        seeds_filename="seeds.dat"):

    def odd_position(column_number, row_number):
        #This funcrtion creates seeds in a row for odd columns

        cell_center = [cell_width * column_number + cell_width/2,  (2 * cell_height * row_number + cell_height/2)]

        left_treshold = cell_center[0] - width_treshold/2
        lower_treshold = cell_center[1] - width_treshold/2

        random_point = [left_treshold + width_treshold * random.random(), lower_treshold + height_treshold * random.random()]
        return(random_point)

    def even_position(column_number, row_number):
        #This funcrtion creates seeds in a row for even columns
        cell_center = [cell_width * column_number, (2* cell_height * row_number + 3*cell_height/2)]

        left_treshold = cell_center[0] - width_treshold/2
        lower_treshold = cell_center[1] - width_treshold/2

        if column_number == 1:
            random_point = [cell_center[0] + width_treshold * random.random()/2, lower_treshold + height_treshold * random.random()]

        elif column_number == column:
            random_point = [left_treshold + width_treshold * random.random()/2, lower_treshold + height_treshold * random.random()]

        else:
            random_point = [left_treshold + width_treshold * random.random(), lower_treshold + height_treshold * random.random()]

        return(random_point)

    domain_x = domain_y * numpy.sqrt(3)/1.5
    column = 2 * row
    
    seeds = numpy.zeros(((2 * column ) * row,2))
    cell_height = domain_y/row/2
    cell_width = domain_x/column 

    width_treshold = cell_width * DoI
    height_treshold = cell_height * DoI

    k = 0
    for i in range(row):
        for j in range(column):
            seeds[2*k] = odd_position(j, i)
            seeds[2*k+1] = even_position(j, i)
            k += 1
    print(seeds)

    gen.write_seeds_to_file(seeds, seeds_filename)

    return seeds
