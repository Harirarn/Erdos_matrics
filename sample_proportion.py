import erdos_finder
import random
import numpy as np

import configparser
config = configparser.ConfigParser()
config.read("config.ini")
N = int(config["Marco-Polo"]["N"])
sample_size = int(config["Marco-Polo"]["sample_size"])

def random_skel(N=N):
    return np.array([[random.randrange(2) for i in range(N)] for j in range(N)])

admis_counter = 0
erdos_counter = 0
admis = []
erdos = []
for i in range(sample_size):
    S = random_skel()
    if erdos_finder.admissible(S):
        admis_counter += 1
        admis.append(S)
        E = erdos_finder.erdosify_from_mat(S)
        if E and E.is_erdos and erdos_finder.same_zeros(S, E.numerator):
            erdos_counter += 1
            erdos.append(E)

print(f"{erdos_counter}:{admis_counter}={erdos_counter/admis_counter}")
        
        
