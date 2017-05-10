module load Langs/Intel/15
module load Langs/Intel/15 MPI/OpenMPI/1.8.6-intel15

mpiexec -n 11 parallel 100 100 10 10 0 10
mpiexec -n 11 parallel 100 100 10 10 1 10
mpiexec -n 11 parallel 100 100 10 10 2 10
mpiexec -n 11 parallel 100 100 10 10 3 10
mpiexec -n 11 parallel 100 100 10 10 4 10


