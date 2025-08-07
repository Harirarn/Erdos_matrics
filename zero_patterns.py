import itertools
import numpy as np
import erdos_finder
import pickle
from time import time

N = 5

identity = np.identity(N, dtype=int)
permutations = [np.array(par) for par in itertools.permutations(identity)]
permutations_lis = [np.array(par) for par in itertools.permutations(np.arange(N))]


complementy = np.ones(N, dtype=int) - identity
pos_values = np.array([[2**(i+N*j) for i in range(N)] for j in range(N)], dtype = int)
ideed_pos_values = np.array([[[2**(i+(N-1)*j - (i>j)), 0][i==j] for i in range(N)] for j in range(N)], dtype = int)
ideed_div = ideed_pos_values + identity
stair_identity = np.diag([2**i for i in range(N*N-N, N*N)])
ideed_pos_complete = ideed_pos_values + stair_identity
id_hash =  np.sum(stair_identity)

def mat_hash(mat):
    return np.sum(mat*pos_values)
##
##def mat_from_num(n):
##    return np.array([[(n//(2**(i+N*j)))%2 for i in range(N)] for j in range(N)])

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
        

def simple_sum(perms):
    A = np.zeros((N, N), int)
    for perm in perms:
        for i in range(len(perm)):
                A[i,perm[i]] += 1
    return A

def simple_average(perms):
    A = simple_sum(perms)
    gcd = np.gcd.reduce(A.flat)
    A //= gcd
    den = int(np.sum(A[0]))
    return erdos_finder.fracMatrix(A, den, None, erdos_finder.linearCombination(np.ones(len(perms))/len(perms), perms))

def same_zeros(A, B):
    return ((A==0) == (B==0)).all()

def admissible(mat):
    for i in range(N):
        if np.sum(mat[i]) == 0 or np.sum(mat[:,i]) == 0:
            return False
    support = erdos_finder.supporting_permutations(mat, permutations_lis)
    if same_zeros(mat, simple_sum(support)):
        return True
    return False

def linearly_independent_combinations(gram, once = True):
    selection = []
    rank = np.linalg.matrix_rank(gram)
    i = 0
    while True:
        if i < len(gram):
            selection.append(i)
            sub_gram_matrix = gram[np.ix_(selection, selection)]
            if not np.linalg.det(sub_gram_matrix) >= 0.1:
                selection.pop()
            i += 1
        else:
            if len(selection) >= rank:
                yield selection
                if once:
                    return
            if not selection:
                return
            i = selection.pop()+1
def basis_from_gram(gram):
    return list(linearly_independent_combinations(gram))[0]

if __name__ == "__main__":

    pickle_file = "erdos" + str(N) + ".pkl"
    try:
        raise
        with open(pickle_file, mode="rb") as f:
            reps, preps, resultants, distinction = pickle.load(f)
    except:
        reps = []
        preps = []
        resultants = {}
        distinction = []

    start = time()
    if not reps:
        partial_filename = f"progress{N}.pkl"
        completed = False
        try:
            total = 2**(N*(N-1))
            print(f"Sieving through {total} boolean matrices.")
            try:
                with open(partial_filename, "rb") as partial_file:
                    initat, reps, reps_if, distinction = pickle.load(partial_file)
                    print(f"Found previous work. Continuing from i = {initat}")
            except FileNotFoundError:
                initat, reps_if = 0, np.ones(total, dtype=bool)
            
                    
            for i in range(initat, total):
                if reps_if[i]:
                    reps.append(i)
                    mat = mat_from_num(i)
                    for eqnum in (equivs:=num_equivalents(mat)):
                        if eqnum >= 0:
                            reps_if[eqnum] =  False
                    distinction.append((len(equivs),mat.sum(),i))
                    
            print(f"{len(reps)} representatives found. {time()-start:0.2f} s passed.")
            completed = True
        except KeyboardInterrupt:
            print("Keyboard innterrupt occured. Saving work")
            with open(partial_filename, "wb") as partial_file:
                pickle.dump((i, reps, reps_if, distinction), partial_file)
            print("Work saved")

    preps = [i for i in reps if admissible(mat_from_num(i))]
    print(f"Pruned it down to {len(preps)} admissible representatives.")
    print(f"{len(reps)} representatives found so far. {time()-start:0.2f} s passed.")

    for n in preps:
        if n >= initat:
            mat = mat_from_num(n)
            support_perms = erdos_finder.supporting_permutations(mat, permutations_lis)

            if support_perms:
                # Calculating the simple average of all matrices in support. In some cases this is Erdos.
                E = simple_average(support_perms)
                if not E.is_erdos: # Running the algorithm if simple average is not good enough
                    indices = basis_from_gram(erdos_finder.gram_matrix(support_perms))
                    selection = [support_perms[j] for j in indices]
                    E = erdos_finder.erdosify(selection)
                resultants[n] = E

    erdoses = {}
    for n, E in resultants.items():
        mat = mat_from_num(n)
        if E and E.is_erdos and same_zeros(mat, E.numerator):
            erdoses[n] = E

    print(f"Found {len(erdoses)} Erdos matrices. {time()-start:0.2f} s passed.")
                    
    distinction.sort()

    with open(pickle_file, mode="wb") as f:
        pickle.dump((reps,preps,resultants,distinction), f)


def dict_minus(first, other=None):
    ret = {}
    if other is None:
        other = first
    for n in first:
        for m in other:
            if(first[n].numerator == other[m].numerator).all():
                break
        else:
            ret[n] = first[n]
    return ret
    
