import numpy as np
import itertools
from dataclasses import dataclass
import functools
from time import time
import pickle


@functools.lru_cache()
def ideed_pos_complete(N):
    ideed_pos_values = np.array([[[2**(i+(N-1)*j - (i>j)), 0][i==j] for i in range(N)] for j in range(N)], dtype = int)
    stair_identity = np.diag([2**i for i in range(N*N-N, N*N)])
    return ideed_pos_values + stair_identity

@functools.lru_cache()
def id_hash(N):
    return 2**(N*N) - 2**(N*N-N) # = np.sum(stair_identity)

@functools.lru_cache()
def num_from_matN(N):
    def num_from_mat(mat, rectified = id_hash(N), idier = ideed_pos_complete(N)):
        return np.sum(mat*idier) - rectified
    return num_from_mat

@functools.lru_cache()
def mat_from_numN(N):
    def mat_from_num(n, rectified = id_hash(N), idier = ideed_pos_complete(N)):
        n += rectified
        return n // idier %2
    return mat_from_num


@functools.lru_cache()
def permutations_np(N):
    return [np.array(par) for par in itertools.permutations(np.identity(N, dtype=int))]

@functools.lru_cache()
def permutations_lis(N):
    return [np.array(par) for par in itertools.permutations(np.arange(N))]

def trace_along(A, perm):
    return sum(A[i, p] for i, p  in enumerate(perm))
def allTrace(A, permlist=None):
    if permlist is None: permlist = permutations_lis(len(A[0]))
    return np.array([trace_along(A, perm) for perm in permlist])
def maxTrace(A, permlist=None):
    """ Evaluates the maximum across trace of matrices with rows permuted."""
    return max(allTrace(A, permlist))

def gram_matrix(A, B = None):
    if B is None: B = A
    return np.array([[sum(a == b) for b in B] for a in A])


@dataclass
class linearCombination:
    coefficients: np.array
    elements: list

    def perm_sum(self):
        N = len(self.elements[0])
        ret_matrix = np.zeros((N, N), dtype = int)
        for coeff, perm in zip(self.coefficients, self.elements):
            for i in range(len(perm)):
                ret_matrix[i,perm[i]] += coeff
        return ret_matrix
        

@dataclass
class fracMatrix:
    numerator: np.array
    denominator: int
    is_erdos: bool = None
    solution: linearCombination = None

    def __post_init__(self):
        if self.is_erdos is None:
            self.verify()

    def matrix(self):
        return self.numerator/self.denominator
    def __str__(self):
        return f"numerator=\n{self.numerator},\ndenominator={self.denominator}, is_erdos={self.is_erdos}"
    def verify(self):
        self.is_erdos = (self.numerator >= 0).all() and bool(np.sum(self.numerator**2) == maxTrace(self.numerator)*self.denominator)

def erdosify(perms):
    """ For a given list of LI permutation matrices, finds an Erdos matrix in their linear span."""
    N = len(perms[0])
    gram = gram_matrix(perms)
    if np.linalg.det(gram) >= 0.1:
        sol = np.linalg.solve(gram, np.ones(len(gram)))
        # make the coefficients into integers
        sol = np.array(np.linalg.det(gram)*sol + 0.01 - 0.02*(sol<0), dtype = np.int64)
        sol //= np.gcd.reduce(sol)
        denominator = sum(sol)
        E = linearCombination(sol, perms)
        numerator = E.perm_sum()
        if (d:=np.gcd.reduce(numerator.flat)) != 1:
            numerator //= d
            denominator //= d
        frob_norm = np.sum(numerator**2)
        is_erdos = (numerator >= 0).all() and (frob_norm == trace_along(numerator, perms[0])*denominator) and (frob_norm == maxTrace(numerator)*denominator)
        return fracMatrix(numerator, int(denominator), bool(is_erdos), E)
    else:
        return False

def erdosify_from_support(support_perms):
    """ Similare to erdosify, but if passed a lin-dependent set, first gets their maximal independent subset."""
    indices = basis_from_gram(gram_matrix(support_perms))
    selection = [support_perms[j] for j in indices]
    E = erdosify(selection)
    return E
def erdosify_from_mat(mat):
    """ Similare to erdosify, but takes a skeleton as its input"""
    support_perms = supporting_permutations(mat)
    return erdosify_from_support(support_perms)
        
def perm_format(perm):
    """ Returns the cycle notation of a permutation. 0 is still the first index."""
    remaining = list(range(len(perm)-1, -1, -1))
    ret = []
    current = []
    while remaining:
        if current:
            next = perm[current[-1]]
            current.append(int(next))
            remaining.remove(next)
        else:
            current = [remaining.pop()]
        if perm[current[-1]] == current[0]:
            if len(current)>1:
                ret.append(current)
            current = []
    return ret

def preservers(perm, other = None):
    N = len(perm)
    if other is None: other = perm
    for x, y in itertools.product(permutations_np(N), repeat=2):
        if (perm[np.ix_(x, y)] == other).all():
            yield (x, y)

# These are now called inner permutations in the paper
def supporting_permutations(matrix, permlist=None):
    """ Returns the list of permutations that for which all entries in matrix is non-zero"""
    if permlist is None: permlist = permutations_lis(len(matrix[0]))
    return [perm for perm in permlist if 0 not in [matrix[i, p] for i, p  in enumerate(perm)]]

def comb(n, r):
    num = den = 1
    r = min(r, n-r)
    for i, j in zip(range(1, r+1), range(n, n-r, -1)):
        num *= j
        den *= i
    return num // den

def perm_format_from_matrix(pnp):
    """ If the passed matrix is a permutation, then returns in cycle notation."""
    N = pnp[0]
    for pl, pn in zip(permutations_lis(N), permutations_np(N)):
        if (pn == pnp).all():
            return perm_format(pl)

def simple_sum(perms):
    N = len(perms[0])
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
    return fracMatrix(A, den, None, linearCombination(np.ones(len(perms))/len(perms), perms))


def same_zeros(A, B):
    return ((A==0) == (B==0)).all()

def admissible(mat):
    """ Returns if the given skeleton is admissible or not. """
    N = len(mat[0])
    for i in range(N):
        if np.sum(mat[i]) == 0 or np.sum(mat[:,i]) == 0:
            return False
    support = supporting_permutations(mat, permutations_lis(N))
    if support and same_zeros(mat, simple_sum(support)):
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

def distinction_stater(n, E):
    return (len(set(E.numerator.flat)),
            list(E.numerator.flat).count(0),
            n)

if __name__ == "__main__":
    import configparser
    config = configparser.ConfigParser()
    config.read("config.ini")
    N=int(config["DEFAULT"]["N"])
    mat_from_num = mat_from_numN(N)
    num_from_mat = num_from_matN(N)
    pickle_file = f"pickles/erdos{N}.pkl"
    try:
        #raise # Uncomment this raise to force fresh execution.
        with open(pickle_file, mode="rb") as f:
            reps, preps, resultants, distinction, erdoses = pickle.load(f)
    except:
        reps = {}
        resultants = {}
        distinction = []

    if not reps:
        # Loading the class representatives from rep file. Run reps_multi_flat to generate this.
        start = time()
        with open(f"pickles/reps{N}.pkl", "rb") as result_file:
            reps = pickle.load(result_file)
        print(f"Found {len(reps)} representatives.")

        # Collecting admissible skeletons
        preps = {i:j for i,j in reps.items() if admissible(mat_from_num(i))}
        print(f"Pruned it down to {len(preps)} admissible representatives. {time()-start:0.2f} s passed.")

        for n in preps:
            mat = mat_from_num(n)
            support_perms = supporting_permutations(mat, permutations_lis(N))

            if support_perms:
                # Calculating the simple average of all permutations in support. In some cases this is Erdos.
                E = simple_average(support_perms)
                if not E.is_erdos: # Running the algorithm when simple average is not good enough
                    E = erdosify_from_support(support_perms)
                resultants[n] = E
                distinction.append(distinction_stater(n, E))

        erdoses = {n:E for n, E in resultants.items() if E and E.is_erdos and same_zeros(mat_from_num(n), E.numerator)}
        distinction.sort()

        print(f"Found {len(erdoses)} Erdos matrices. {time()-start:0.2f} s passed.")
                    

        with open(pickle_file, mode="wb") as f:
            pickle.dump((reps,preps,resultants,distinction,erdoses), f)

    # Creating a .txt file with erdos matrices.
    with open(f"erdos{N}x{N}.txt", "w") as f:
        for E in erdoses.values():
            f.write(str(E) + "\n\n")
