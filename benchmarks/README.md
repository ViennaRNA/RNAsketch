Hoehner Multistate example Benchmark
====================================

These benchmarks runs Hoehner tristable example and plot how many iteration of several options we need
to reach a minimal difference between eos(state) and the mfe
It also plots how often we get the mfe being one of the input states and the score;

 * random_optimization.py very simple, just sample sequences and remember best solution. sampling can be done completely random
 or with mutation and therefore along a path through the solution space
 * exit_duration.py benchmark for different number of neutral trials before exiting the optimization
 * jump_duration.py benchmark for different number of random jumps before the actual optimization


Results plotting:
-----------------


Random Sequence -> Subopt -> Design Benchmark
=============================================

Here we get a random sequence, run a RNAsubopt, take some suboptimal structures
and the mfe structure and try to get the initial Sequence again.

