# Task-Pool-management
Dynamic and Static Task Pool management
The purpose of the exercise is to implement a simple application with Dynamic and Static Task Pool management.
Using MPI with C.

To run the program: (in linux)
1. Open terminal in the file solution
2. run: mpicc HW.c -o HW -lm
        mpiexec -np <numberOfProcesses> ./HW

For n slaves please run n+1 processes (+1 is for the master).
