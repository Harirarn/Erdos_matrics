from multiprocessing import Process, Queue, Pool
import queue
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

def num_equivalents_worker(job_queue, return_queue):
    while True:
        n = job_queue.get()
        if n is None: break
        ret = num_equivalents(mat_from_num(n))
        return_queue.put((n, ret))

def reps_if_manager(job_queue, return_queue, pool_size = 4, batch_size = 1):
    total = 2**(N*(N-1))
    reps_if = np.ones(total, dtype=bool)
    reps = []
    start = time()
    allocated = 0
    waiting_time = 0
    wasted = 0
    i = 0
    increment = (2**N-1)*(2**(N-1)-1)
    while time()-start < 250:
        # Post a job
        while allocated <= batch_size*pool_size:
            if i >= total:
                allocated += 1
            elif reps_if[i]:
                job_queue.put(i)
                allocated += 1
            if i < total:
                i = (i+increment)%total
            if i == 0:
                i += total
        # Process a returned job
        try:
            wait = time()
            #print(job_queue.qsize(), return_queue.qsize())
            n, ret = return_queue.get(timeout=5)
            waiting_time += time()-wait
            allocated -= 1
            if reps_if[n]:
                reps.append(n)
                for eqnum in ret:
                    if eqnum >= 0:
                        reps_if[eqnum] = False
            else:
                wasted += 1
        except queue.Empty:
            break
    print(i, waiting_time, wasted)
    print(time()-start)
    return (reps, (len(reps_if)-reps_if.sum())/len(reps_if))            
        

if __name__ == "__main__":
    job_queue = Queue()
    return_queue = Queue()
    cores = 6
    workers = []
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
    for p in range(cores):
        w = workers.pop()
        w.kill()
        #w.join()
        #w.close()
