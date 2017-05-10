module load Langs/Intel/15
module load Langs/Intel/15 MPI/OpenMPI/1.8.6-intel15

echo "format: serial/approx/parallel time/fitness" echo "10 x 5"
mpiexec -n 11 parallel 10 10 10 5 5 10
mpiexec -n 11 parallel 100 100 10 5 5 10
mpiexec -n 11 parallel 1000 100 10 5 5 10
mpiexec -n 11 parallel 2000 200 10 5 5 10
mpiexec -n 11 parallel 5000 500 10 5 5 10
mpiexec -n 11 parallel 10000 1000 10 5 5 10


