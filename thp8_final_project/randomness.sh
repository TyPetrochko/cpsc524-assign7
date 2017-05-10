module load Langs/Intel/15
module load Langs/Intel/15 MPI/OpenMPI/1.8.6-intel15
echo "format: serial_time   serial_fitness    approx_time   approx_fitness  parallel_time    parallel_fitness"

mpiexec -n 11 parallel 100 100 10 10 0 10
mpiexec -n 11 parallel 100 100 10 10 1 10
mpiexec -n 11 parallel 100 100 10 10 2 10
mpiexec -n 11 parallel 100 100 10 10 3 10
mpiexec -n 11 parallel 100 100 10 10 4 10
mpiexec -n 11 parallel 100 100 10 10 5 10
mpiexec -n 11 parallel 100 100 10 10 6 10
mpiexec -n 11 parallel 100 100 10 10 7 10


