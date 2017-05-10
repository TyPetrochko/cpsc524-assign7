module load Langs/Intel/15
module load Langs/Intel/15 MPI/OpenMPI/1.8.6-intel15

echo "format: serial/approx/parallel time/fitness"
echo "10 x 5"
mpiexec -n 11 parallel 1000 100 10 5 4 10
mpiexec -n 11 parallel 1000 1000 10 5 4 10
mpiexec -n 11 parallel 1000 10000 10 5 4 10

echo "10 x 10"
mpiexec -n 11 parallel 1000 100 10 10 4 10
mpiexec -n 11 parallel 1000 1000 10 10 4 10
mpiexec -n 11 parallel 1000 10000 10 10 4 10

echo "10 x 40"
mpiexec -n 11 parallel 1000 100 10 40 4 10
mpiexec -n 11 parallel 1000 1000 10 40 4 10
mpiexec -n 11 parallel 1000 10000 10 40 4 10

