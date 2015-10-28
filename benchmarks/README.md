Hoehner Multistate example Benchmark
====================================

This benchmark runs Hoehner tristable example and plots how many optimization iterations we need for 
global or local optimization to reach a minimal difference between eos(state) and the mfe
It also plots how often we get the mfe being one of the input states.

optimization_duration.py
grid_start_optimization_global.sh
grid_start_optimization_local.sh
results go to: optimization_duration_output

Results plotting:
-----------------


Random Sequence -> Subopt -> Design Benchmark
=============================================

Here we get a random sequence, run a RNAsubopt, take some suboptimal structures
and the mfe structure and try to get the initial Sequence again.



