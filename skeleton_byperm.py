import itertools
import numpy as np
import pickle
from time import time

N = 6

identity = np.identity(N, dtype=int)
permutations_lis = [np.array(par) for par in itertools.permutations(np.arange(N))]


ideed_pos_values = np.array([[[2**(i+(N-1)*j - (i>j)), 0][i==j] for i in range(N)] for j in range(N)], dtype = int)
stair_identity = np.diag([2**i for i in range(N*N-N, N*N)])
ideed_pos_complete = ideed_pos_values + stair_identity
## The above matrix will be of the form
## 2** [N(N-1) 0 1 2 3 ...       N-2]
##     [N-1  N(N-1)+1 N N+1 ... 2N-3]
##     [ ...                        ]
##     [ ...                  N**2-1]
id_hash =  np.sum(stair_identity)

def num_from_mat(mat, rectified = id_hash, idier = ideed_pos_complete):
    return np.sum(mat*idier) - rectified

def mat_from_num(n, rectified = id_hash, idier = ideed_pos_complete):
    n += rectified
    return n // idier %2

def perm_equivalents(mat):
    for L, R in itertools.product(permutations_lis, repeat=2):
        yield mat[np.ix_(L, R)]

ideed_equivalents = list(perm_equivalents(ideed_pos_complete))

def num_equivalents(mat, rectified = id_hash):
    eqnums = {num_from_mat(mat, rectified, idier) for idier in ideed_equivalents}
    if num_from_mat(mat.T) not in eqnums:
        eqnums.update(num_from_mat(mat.T, rectified, idier) for idier in ideed_equivalents)
    return eqnums

if __name__ == "__main__":
    partial_filename = f"progress{N}.pkl"
    try:
        total = 2**(N*(N-1))
        print(f"Sieving through {total} boolean matrices.")
        try:
            with open(partial_filename, "rb") as partial_file:
                initat, reps, reps_if, distinction = pickle.load(partial_file)
                print(f"Found previous work. Continuing from i = {initat}")
        except FileNotFoundError:
            initat, reps, reps_if, distinction = 0, [], np.ones(total, dtype=bool), []
        
                
        for i in range(initat, total):
            if reps_if[i]:
                reps.append(i)
                mat = mat_from_num(i)
                for eqnum in (equivs:=num_equivalents(mat)):
                    if eqnum >= 0:
                        reps_if[eqnum] =  False
                distinction.append((len(equivs),mat.sum(),i))
                
        print(f"{len(reps)} representatives found. {time()-start:0.2f} s passed.")
        with open(f"reps{N}.pkl", "wb") as result_file:
            pickle.dump((reps, distinction), result_file)
    except KeyboardInterrupt:
        print("Keyboard innterrupt occured. Saving work")
        with open(partial_filename, "wb") as partial_file:
            pickle.dump((i, reps, reps_if, distinction), partial_file)
        print("Work saved. {time()-start:0.2f} s passed.")
