This code is concerned with finding Erdos matrices. Check https://arxiv.org/abs/2512.04766 for relavant theory behind this.

# Relavant files

The .py files require `python >=3.10` and `numpy` installed.
## config.ini
- N : The other listed files will work with NxN matrices. Filenames containing {N} below are this value.
- cores : Files using multiprocessing will use this many cores.
- max_time : The amount of time in seconds for which reps_multi_flat.py will run per session. Set this to 0 to remove this limit.

## reps_multi_flat.py
```sh
python3 reps_multi_flat.py
```
Outputs: pickles/reps{N}.pkl, pickles/progress{N}.pkl

This executes step 1 of the algorithm, i.e. finds a list of representative skeletons for the equivalence classes of skeletons under permutation of rows and columns and under transposition.
This outputs as a dict in reps{N}.pkl pickle file. The dict follows 'index of representative skeleton S':('size of its equivalence class', 'if symmetrizable by permutation of rows and columns').
This program is also capable of storing and resuming from partial progress of the task. It runs for max_time (from config.ini) seconds per session. It should be safe to send a single KeyboardInterrupt during execution and progress will be saved. The progress is stored in progress{N}.pkl. Warning: progress{N}.pkl requires 1 GB disk space when N=6.
Dont run this file for N>6.

## erdos_finder.py
```sh
python3 erdos_finder.py
```
Requires: pickles/reps{N}.pkl

Outputs: pickles/erdos{N}.pkl, erdos{N}x{N}.txt

This program executes rest of the steps in the algorithm. The output is stored in erdos{N}.pkl as a tuple (reps, preps, resultants, erdoses, distinction).
- reps: this is the same dict that was in reps{N}.pkl.
- preps: this is a sub-dict of reps with indices only when the skeleton is admissible.
- resultants: this is a dict with keys from preps (indices of admissible skeletons) and values now representing the final result of the algorithm (whether erdos or not).
- erdoses: this is a sub-dict of resultants with indices only when the matrix is erdos.
- distinction: this is a list of tuples (one for each entry in resultants) as ('number of distinct entries', 'number of zeros', 'index'). This also gets sorted in increasing order.
If you wish to unpickle this file, add the following line to ensure relevant dataclass are defined in that scope:
```python
from erdos_finder import fracMatrix, linearCombination
```

The list of erdos matrices are also outputted into the text file erdos{N}x{N}.txt

## erdos_failure_stats.py
```sh
python3 erdos_failure_stats.py
python3 erdos_failure_stats.py >stats{N}.txt
```
Requires: pickles/erdos{N}.pkl

This program evaluates and prints a bunch of relavant stats. They have also been recorded in stats{N}.txt


Other .py files are old ones and may be ignored for now.
