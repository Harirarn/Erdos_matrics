from multiprocessing import Process, Queue, Pool
import queue
import itertools
import numpy as np
import pickle
from time import time

import configparser
config = configparser.ConfigParser()
config.read("config.ini")
N=int(config["DEFAULT"]["N"])
cores=int(config["DEFAULT"]["cores"])

identity = np.identity(N, dtype=int)
permutations_lis = [np.array(par) for par in itertools.permutations(np.arange(N))]


ideed_pos_values = np.array([[[2**(i+(N-1)*j - (i>j)), 0][i==j] for i in range(N)] for j in range(N)], dtype = int)
stair_identity = np.diag([2**i for i in range(N*N-N, N*N)])
ideed_pos_complete = ideed_pos_values + stair_identity
# The above matrix is used for place values.
# Log_2(ideed_pos_complete) is
# [[ N**2-N 0 1 2 ... N-2]
#  [ N-1 N**2-N+1 N N+1 ... 2*N-3]
# ...
#  [ N**2-2*N+1 ... N**2-N-1 N**2-1]]
# The positions 0, ..., N**2-N-1 are used for the offdiagonal entries and the
# remaining N**2-N, ... N**2-1 are used for the diagonal.


ideed_pos_complete_trans = ideed_pos_complete.T
ideed_pos_complete_flat = ideed_pos_complete.flatten()
ideed_pos_complete_trans_flat = ideed_pos_complete_trans.flatten()
# The matrices are flattened in this file as it speeds up the numpy computations.
id_hash =  np.sum(stair_identity)

# The following two functions translate back and forth between the binary string and matrix.
def num_from_mat(mat, idier = ideed_pos_complete_flat, rectified = id_hash):
    return np.dot(mat, idier) - rectified

def mat_from_num(n, idier = ideed_pos_complete_flat, rectified = id_hash):
    n += rectified
    return n // idier %2

def perm_equivalents(mat):
    for L, R in itertools.product(permutations_lis, repeat=2):
        yield mat[np.ix_(L, R)]

ideed_equivalents = [idier.flatten() for idier in perm_equivalents(ideed_pos_complete)]
ideed_equivalents_trans = [idier.flatten() for idier in perm_equivalents(ideed_pos_complete_trans)]

def num_equivalents(n, rectified = id_hash):
    """ Returns a set containing the id of all equivalent matrices, and if symmetrizable by permutations."""
    mat_flat = mat_from_num(n)
    eqnums = {num_from_mat(mat_flat, idier, rectified) for idier in ideed_equivalents}
    if (trans_asym:=(num_from_mat(mat_flat, ideed_pos_complete_trans_flat, rectified) not in eqnums )):
        eqnums.update(num_from_mat(mat_flat, idier, rectified) for idier in ideed_equivalents_trans)
    return eqnums, not trans_asym

def num_equivalents_worker(job_queue, return_queue):
    while True:
        n = job_queue.get()
        if n is None: break
        ret = num_equivalents(n)
        return_queue.put((n, ret))

def reps_if_manager(job_queue, return_queue, pool_size = 4, batch_size = 1):
    total = 2**(N*(N-1))
    try: # We check if we had some partial progress done in a previous session.
        #raise FileNotFoundError # Uncomment this line to force fresh execution. Alternatively delete the progress file.
        with open(f"progress{N}.pkl", "rb") as partial_file:
            i, ni, allocated, reps, reps_if = pickle.load(partial_file)
        for n in allocated:
            job_queue.put(n)
    except FileNotFoundError:
        i = 0
        ni = 0
        allocated = []
        reps = {}
        reps_if = np.ones(total, dtype=bool)
    
    start = time()
    waiting_time = 0 
    wasted = 0
    # Neighbouring matrices are often in the same equivalence class.
    # By incrementing by a bigger number we ensure less wasted work.
    # The quantity (2**N-1)*(2**(N-1)-1) is a divisor of 2**(N(N-1))-1. This ensures that we go through the whole list.
    increment = (2**N-1)*(2**(N-1)-1) 
    max_allocations = batch_size*pool_size+1
    try:
        while time()-start < 28800: # If you want to keep the computer running for longer per session, increase the number or remove this check.
            # Post a job
            while len(allocated) < max_allocations:
                if i >= total:
                    job_queue.put(None)
                    allocated.append(0)
                elif reps_if[ni]:
                    job_queue.put(ni)
                    allocated.append(ni)
                if i < total:
                    ni = (ni+increment)%total
                    i+=1
            # Process a returned job
            try:
                wait = time()
                #print(job_queue.qsize(), return_queue.qsize()) #Debug
                n, ret = return_queue.get(timeout=20)
                waiting_time += time()-wait
                if reps_if[n]:
                    reps[n] = (len(ret[0]), ret[1])
                    for eqnum in ret[0]:
                        if eqnum >= 0:
                            reps_if[eqnum] = False
                else:
                    wasted += 1
                allocated.remove(n)
                if i == total and sum(allocated)==0 and allocated:
                    break
            except queue.Empty:
                waiting_time += time()-wait
                break
    finally:
        print("saving progress.")
        with open(f"progress{N}.pkl", "wb") as partial_file:
            pickle.dump((i, ni, allocated, reps, reps_if), partial_file)
        print("saving complete")
            
        print(f"index reached:{i}\n manager_waiting_time:{waiting_time}\n wasted_job_allocations:{wasted}")
        print(f"{time()-start}s passed.")
    return (reps, (len(reps_if)-reps_if.sum())/len(reps_if))            
        

if __name__ == "__main__":
    job_queue = Queue()
    return_queue = Queue()
    workers = []
    try:
        for i in range(cores):
            workers.append(Process(target = num_equivalents_worker, args = (job_queue, return_queue)))
            workers[-1].start()
        result = reps_if_manager(job_queue, return_queue, cores)
        reps, done = result
        print(len(reps), done)
        with open(f"reps{N}.pkl", "wb") as result_file:
            pickle.dump(reps, result_file)
        for p in range(cores+1):
            job_queue.put(None)
    finally:
        while workers:
            w = workers.pop()
            w.kill()
        #w.join()
        #w.close()
