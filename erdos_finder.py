import numpy as np
import itertools
from dataclasses import dataclass

N = 3

def gram(parms):
    parms

identity_np = np.identity(N, dtype=int)
identity_lis = np.arange(N)

permutations_np = [np.array(par) for par in itertools.permutations(identity_np)]
permutations_lis = [np.array(par) for par in itertools.permutations(identity_lis)]

other_permutations = permutations_np[1:]

def maxTrace(A, permlist=permutations_lis):
    """ Evaluates the maximum across trace of matrices with rows permuted."""
    return max(sum(A[i, perm[i]] for i in range(N)) for perm in permlist)


def splTrace(A):
    """ A heuristic function which helps in avoiding calculating maxTrace. Not much saving for N <= 4"""
    u = len(A)
    B = A.copy()
    ret = 0
    for i in range(u): #Iteratively finds the max entry and removes its column and row. Adds all these entries up.
        max_entry = max(B.flat)
        ret += max_entry
        indices = np.where(B==max_entry)
        maxi, maxj = indices[0][0], indices[1][0]
        for j in range(u):
            B[maxi, j] = B[j, maxj] = 0
    return ret

max_degree = (N-1)**2+1
max_degree = 2

def gramMatrix(A, B = None):
    if B is None: B = A
    return np.array([[sum(a == b) for b in B] for a in A])


@dataclass
class fracMatrix:
    numerator: any
    denominator: int
    is_erdos: bool

    def matrix(self):
        return self.numerator/self.denominator

def erdosify(perms):
    gram = gramMatrix(perms)
    if np.linalg.det(gram) != 0:
        sol = np.linalg.solve(gram, np.ones(len(gram)))
        # make the coefficients into integers
        sol = np.array(np.linalg.det(gram)*sol + 0.01, dtype = np.int64)
        sol //= np.gcd.reduce(sol)
        denominator = sum(sol)
        ret_matrix = np.zeros((N, N), dtype = int)
        for coeff, perm in zip(sol, perms):
            for i in range(N):
                ret_matrix[i,perm[i]] += coeff

        is_erdos = maxTrace(ret_matrix)*denominator == np.sum(ret_matrix**2)
        return fracMatrix(ret_matrix, int(denominator), bool(is_erdos))
        


##for k in range(max_degree):
##    for selection in itertools.combinations(permutations_lis[1:], k):
##        perms = list(selection)
##        perms.append(identity_lis)
##        erdosify(perms)

if __name__ == "__main__":
    E = erdosify(permutations_lis[:2])
    print(E)
