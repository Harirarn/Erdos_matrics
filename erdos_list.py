from dataclasses import dataclass
from erdos_finder import *

@dataclass
class fracMatrix:
    numerator: any
    denominator: int
    is_erdos: bool

    def matrix(self):
        return self.numerator/self.denominator

erdos4 = [
    fracMatrix(np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]], dtype = int), 1, True),
    fracMatrix(np.array([[1,0,0,1],[0,2,0,0],[0,0,2,0],[1,0,0,1]], dtype = int), 2, True),
    fracMatrix(np.array([[1,0,0,1],[0,2,0,0],[1,0,1,0],[0,0,1,1]], dtype = int), 2, True),
    fracMatrix(np.array([[1,0,2,2],[0,5,0,0],[2,0,3,0],[2,0,0,3]], dtype = int), 5, True),
    fracMatrix(np.array([[1,0,1,2],[0,4,0,0],[2,0,2,0],[1,0,1,2]], dtype = int), 4, True),
    fracMatrix(np.array([[1,0,0,1],[0,1,1,0],[1,0,1,0],[0,1,0,1]], dtype = int), 2, True),
    fracMatrix(np.array([[1,0,0,1],[0,1,1,0],[0,1,1,0],[1,0,0,1]], dtype = int), 2, True),
    fracMatrix(np.array([[1,0,1,1],[0,3,0,0],[1,0,1,1],[1,0,1,1]], dtype = int), 3, True),
    fracMatrix(np.array([[1,0,0,2],[1,2,0,0],[1,0,2,0],[0,1,1,1]], dtype = int), 3, True),
    fracMatrix(np.array([[4,0,9,10],[10,13,0,0],[9,0,14,0],[0,10,0,13]], dtype = int), 23, True),
    fracMatrix(np.array([[1,1,3,3],[4,4,0,0],[0,3,5,0],[3,0,0,5]], dtype = int), 8, True),
    fracMatrix(np.array([[1,0,3,4],[3,5,0,0],[3,0,5,0],[1,3,0,4]], dtype = int), 8, True),
    fracMatrix(np.array([[2,0,3,6],[4,7,0,0],[5,0,6,0],[0,4,2,5]], dtype = int), 11, True),
    fracMatrix(np.array([[2,7,15,19],[19,24,0,0],[15,0,28,0],[7,12,0,24]], dtype = int), 43, True),
    fracMatrix(np.array([[1,0,1,2],[2,2,0,0],[0,2,2,0],[1,0,1,2]], dtype = int), 4, True),
    fracMatrix(np.array([[1,0,1,2],[1,2,1,0],[2,0,2,0],[0,2,0,2]], dtype = int), 4, True),
    fracMatrix(np.array([[3,6,6,14],[13,16,0,0],[13,0,16,0],[0,7,7,15]], dtype = int), 29, True),
    fracMatrix(np.array([[3,0,13,13],[14,15,0,0],[6,7,16,0],[6,7,0,16]], dtype = int), 29, True),
    fracMatrix(np.array([[2,0,2,3],[0,4,3,0],[2,3,2,0],[3,0,0,4]], dtype = int), 7, True),
    fracMatrix(np.array([[1,3,3,7],[6,8,0,0],[6,0,8,0],[1,3,3,7]], dtype = int), 14, True),
    fracMatrix(np.array([[1,1,6,6],[7,7,0,0],[3,3,8,0],[3,3,0,8]], dtype = int), 14, True),
    fracMatrix(np.array([[4,4,5,9],[11,11,0,0],[7,7,8,0],[0,0,9,13]], dtype = int), 22, True),
    fracMatrix(np.array([[4,0,7,11],[9,13,0,0],[5,9,8,0],[4,0,7,11]], dtype = int), 22, True),
    fracMatrix(np.array([[2,0,3,5],[5,5,0,0],[3,3,4,0],[0,2,3,5]], dtype = int), 10, True),
    fracMatrix(np.array([[2,3,5,9],[9,10,0,0],[5,6,8,0],[3,0,6,10]], dtype = int), 19, True),
    fracMatrix(np.array([[2,2,5,9],[9,9,0,0],[5,5,8,0],[2,2,5,9]], dtype = int), 18, True),
    fracMatrix(np.array([[5,0,5,7],[6,5,6,0],[6,5,6,0],[0,7,0,10]], dtype = int), 17, True),
    fracMatrix(np.array([[1,1,1,1],[0,2,0,2],[2,0,2,0],[1,1,1,1]], dtype = int), 4, True),
    fracMatrix(np.array([[1,0,1,2],[1,2,1,0],[1,2,1,0],[1,0,1,2]], dtype = int), 4, True),
    fracMatrix(np.array([[1,2,2,3],[2,3,3,0],[2,3,3,0],[3,0,0,5]], dtype = int), 8, True),
    fracMatrix(np.array([[3,3,4,4],[7,7,0,0],[0,4,5,5],[4,0,5,5]], dtype = int), 14, True),
    fracMatrix(np.array([[3,0,4,7],[4,5,5,0],[4,5,5,0],[3,4,0,7]], dtype = int), 14, True),
    fracMatrix(np.array([[2,3,4,4],[6,7,0,0],[2,3,4,4],[3,0,5,5]], dtype = int), 13, True),
    fracMatrix(np.array([[2,2,3,6],[4,4,5,0],[4,4,5,0],[3,3,0,7]], dtype = int), 13, True),
    fracMatrix(np.array([[1,0,1,1],[0,1,1,1],[1,1,1,0],[1,1,0,1]], dtype = int), 3, True),
    fracMatrix(np.array([[1,1,2,2],[3,3,0,0],[1,1,2,2],[1,1,2,2]], dtype = int), 6, True),
    fracMatrix(np.array([[1,1,1,3],[2,2,2,0],[2,2,2,0],[1,1,1,3]], dtype = int), 6, True),
    fracMatrix(np.array([[2,3,3,3],[3,4,0,4],[3,4,4,0],[3,0,4,4]], dtype = int), 11, True),
    fracMatrix(np.array([[2,2,3,3],[2,2,3,3],[3,3,4,0],[3,3,0,4]], dtype = int), 10, True),
    fracMatrix(np.array([[2,2,2,3],[2,2,2,3],[3,3,3,0],[2,2,2,3]], dtype = int), 9, True),
    fracMatrix(np.array([[1,1,1,1],[1,1,1,1],[1,1,1,1],[1,1,1,1]], dtype = int), 4, True),
]

def perm_product(matrix, perm):
    counter = 1
    for row, index in zip(matrix, perm):
        counter *= row[index]
    return counter

pos_values = np.array([[2**(i+N*j) for i in range(N)] for j in range(N)])
def num_from_mat(mat):
    return np.sum(mat*pos_values)

def mat_from_num(n):
    return np.array([[(n//(2**(i+N*j)))%2 for i in range(N)] for j in range(N)])

def min_rep_num(mat):
    mat = np.array(np.array(mat, dtype = bool), dtype = int)
    tmat = mat.T
    return min(min(num_from_mat(L@mat@R), num_from_mat(L@tmat@R)) for L, R in itertools.product(permutations_np, repeat=2))


if __name__ == "__main__":
    test = 2
    if test == 1:
        for F in erdos4:
            at = allTrace(F.numerator)
            maxmzing_perms = [i for i,j in zip(permutations_lis, (at == max(at))) if j]
            support_perms = [p for p in permutations_lis if perm_product(F.numerator, p)>0]
            if len(maxmzing_perms) != len(support_perms):
                print(F)
                print(maxmzing_perms)
                print(support_perms)
    if test == 2:
        indices = [min_rep_num(E.numerator) for E in erdos4]
