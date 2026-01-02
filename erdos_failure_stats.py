import numpy as np
from time import time
import pickle
import erdos_finder
from erdos_finder import fracMatrix, linearCombination


def has_negative(E):
    return (E.numerator<0).any()

def outer_maxtrace(E):
    m = E.numerator
    return erdos_finder.maxTrace(m)*E.denominator > np.sum(m*m)

def smaller_skel(n, E):
    return not erdos_finder.same_zeros(mat_from_num(n), E.numerator)

def dims(E):
    support = erdos_finder.supporting_permutations(E.numerator)
    gram = erdos_finder.gram_matrix(support)
    rank = np.linalg.matrix_rank(gram)
    return rank, len(support)

if __name__ == "__main__":
    import configparser
    config = configparser.ConfigParser()
    config.read("config.ini")
    N=int(config["DEFAULT"]["N"])
    mat_from_num = erdos_finder.mat_from_numN(N)

    with open(f"pickles/erdos{N}.pkl", mode="rb") as f:
        reps, preps, resultants, distinction, erdoses = pickle.load(f)

    print(f"Out of {len(preps)} classes of admissible skeletons, {len(erdoses)} are Erdos.")
    print(f"Ratio: {len(erdoses)/len(preps):.4f}")
    print(f"The reasons for failure of Erdosness of the remaining {len(preps)-len(erdoses)} are:")
    
    failure_state = {n:(has_negative(E), outer_maxtrace(E), smaller_skel(n, E))
                     for n, E in resultants.items()}
    print("Negative entries:", len([s for s in failure_state.values() if s[0]]))
    print("Greater outer trace:", len([s for s in failure_state.values() if s[1]]))
    print("Smaller skeleton:", len([s for s in failure_state.values() if s[2]]))
    print("Negative entries and greater outer trace:", len([s for s in failure_state.values() if s[0] and s[1]]))
    print("Negative entries and smaller skeleton:", len([s for s in failure_state.values() if s[0] and s[2]]))
    print("Greater outer trace and smaller skeleton:", len([s for s in failure_state.values() if s[1] and s[2]]))
    print("All 3 reasons:", len([s for s in failure_state.values() if s[0] and s[1] and s[2]]))

    dims_calc = {n:dims(E) for n, E in resultants.items()}
    simplices = [n for n, d in dims_calc.items() if d[0] == d[1]]
    simplex_erdoses = [n for n in simplices if n in erdoses]
    print(f"\nThere are {len(simplices)} simplices out of which {len(simplex_erdoses)} are erdos")

    print("\nDiscounting the equivalence under transposition...")
    print("Classes of admissible skeletons are:", sum(2-reps[n][1] for n in preps))
    print("Classes of Erdos matrices are:", sum(2-reps[n][1] for n in erdoses))

    print("\nDiscounting all equivalences...")
    print("Number of admissible skeletons are:", admis:=sum(reps[n][0] for n in preps))
    print("Number of Erdos matrices are:", erdos:=sum(reps[n][0] for n in erdoses))
    print(f"Ratio: {erdos/admis:.4f}")

    print("\nErdos matrix with largest denominator:")
    maxdenn = max((E.denominator, n) for n, E in erdoses.items())[1]
    maxdenE = erdoses[maxdenn]
    print(maxdenE)
    
