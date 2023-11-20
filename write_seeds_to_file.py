################################################################################
###                                                                          ###
### Created by Mahdi Manoochehrtayebi, 2021-2022                             ###
### Supervisors: Dr. Aline Bel-Brunon, Dr. Martin Genet                      ###
###                                                                          ###
### École Polytechnique, Palaiseau, France                                   ###
###                                                                          ###
################################################################################

import pickle

################################################################################

def write_seeds_to_file(seeds, filename):
    file = open(filename, "wb")
    pickle.dump(seeds, file)
    file.close()
