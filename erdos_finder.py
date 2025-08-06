import numpy as np
import itertools
from dataclasses import dataclass
import functools

N = 5

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

def supporting_permutations(matrix, permlist=None):
    """ Returns the list of permutations that for which all entries in matrix is non-zero"""
    if permlist is None: permlist = permutations_lis(len(matrix[0]))
    return [perm for perm in permlist if 0 not in [matrix[i, p] for i, p  in enumerate(perm)]]

def comb(n, r):
    num = den = 1
    for i, j in zip(range(1, r+1), range(n, n-r, -1)):
        num *= j
        den *= i
    return num // den

def perm_format_from_matrix(pnp):
    N = pnp[0]
    for pl, pn in zip(permutations_lis(N), permutations_np(N)):
        if (pn == pnp).all():
            return perm_format(pl)

if __name__ == "__main__":
    if N == 4:
        F = fracMatrix(np.array([[3,3,4,4],[7,7,0,0],[0,4,5,5],[4,0,5,5]]),14,True)
        G = fracMatrix(np.array([[3,0,4,7],[4,5,5,0],[4,5,5,0],[3,4,0,7]]),14,True)
        H = fracMatrix(np.array([[2,2,2,3],[2,2,2,3],[3,3,3,0],[2,2,2,3]]),9,True)
        at = allTrace(F.numerator)
        f = [i for i,j in zip(permutations_lis(N), (at == max(at))) if j]
        sup = supporting_permutations(F)
        ex = 2
        for selection in itertools.combinations(f[:ex]+f[ex:], 5):
            perms = list(selection)
            perms.append(f[ex])
            E = erdosify(perms)
            if np.all(E.numerator == F.numerator):
            #if E and E.is_erdos:
                print([perm_format(p) for p in perms])
                print(E)
    elif N == 5:
        T = np.array([[0,0,0,1,1], [0,0,1,1,1], [0,1,1,0,1], [1,1,1,1,0], [1,1,1,0,1]])
        permset = [np.array([0,2,4,3,1]),
                   np.array([1,4,2,3,0]),
                   np.array([2,1,4,3,0]),
                   np.array([3,4,2,0,1]),
                   np.array([4,3,0,2,1]),
                   np.array([4,3,1,2,0]),
            ]
        E = erdosify(permset)
    elif N == 6:
        T = np.array([[1,1,0,0,0,0],
                      [1,1,1,0,0,0],
                      [1,1,1,1,0,0],
                      [1,1,1,1,1,0],
                      [1,0,1,1,0,1],
                      [1,1,1,1,1,1]])
    if N >= 5:
        from time import time
        start_time = time()
        support = supporting_permutations(T)
        print(r:=np.linalg.matrix_rank(gram_matrix(support)))
        c = len(support)
        r -= 0
        total = comb(c, r)
        uniques = []
        print(c)
        for counter, selection in enumerate(itertools.combinations(support, r)):
            E = erdosify(selection)
            if E and E.is_erdos:
                if all(not (E.numerator==U).all() for U in uniques):
                    uniques.append(E.numerator)
                    print([perm_format(p) for p in selection])
                    print(E)
                    print(counter, '/', total)
            if counter%100000 == 0:
                print(counter, '/', total, f"and {time()-start_time:.2f}s elapsed")

