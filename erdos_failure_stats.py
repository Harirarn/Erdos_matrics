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

if __name__ == "__main__":
    import configparser
    config = configparser.ConfigParser()
    config.read("config.ini")
    N=int(config["DEFAULT"]["N"])
    mat_from_num = erdos_finder.mat_from_numN(N)

    with open(f"erdos{N}.pkl", mode="rb") as f:
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

    print("\nDiscounting the equivalence under transposition...")
    print("Classes of admissible skeletons are:", sum(2-p[1] for p in preps.values()))
    print("Classes of Erdos matrices are:", sum(2-preps[n][1] for n in erdoses))

    print("\nDiscounting all equivalences...")
    print("Number of admissible skeletons are:", admis:=sum(p[0] for p in preps.values()))
    print("Number of Erdos matrices are:", erdos:=sum(preps[n][0] for n in erdoses))
    print(f"Ratio: {erdos/admis:.4f}")

    print("\nErdos matrix with largest denominator:")
    maxdenn = max((E.denominator, n) for n,E in erdoses.items())[1]
    maxdenE = erdoses[maxdenn]
    print(maxdenE)
    
