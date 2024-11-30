The program used for the lab report is d1.cpp

In order to reproduce the results:
- Connect to the unitn HPC cluster
- Copy d1.cpp and PBS_d1.pbs in a folder.
- Open PBS_d1.pbs and change the working directory.
- Use qsub PBS_d1.pbs to start the job. 

A d1.o file will be generated containing the printed output, along with details about
the compiler version and the computing architecture of the node.


In the d1.o file, if you also want to visualize the execution times of the functions for 
each iteration, instead of just the computed average, open d1.cpp and remove comments 
from the calls to these 2 functions: print_details_serial and print_details_parallel.

The d1_Output.o file in this folder contains the results presented in the lab report.