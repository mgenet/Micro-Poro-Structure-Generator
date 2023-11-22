################################################################################
###                                                                          ###
### Created by Mahdi Manoochehrtayebi, 2021-2022                             ###
### Supervisors: Dr. Aline Bel-Brunon, Dr. Martin Genet                      ###
###                                                                          ###
### École Polytechnique, Palaiseau, France                                   ###
###                                                                          ###
################################################################################

import numpy

import micro_poro_structure_generator as gen

################################################################################

def generate_seeds_random(
        nb_seeds,
        seeds_filename="seeds.dat"):

    seeds = numpy.random.rand(nb_seeds, 2)

    gen.write_seeds_to_file(seeds, seeds_filename)

    return seeds
