module load Langs/Intel/15
module load Langs/Intel/15 MPI/OpenMPI/1.8.6-intel15

echo "format: serial/approx/parallel time/fitness"
echo "10 x 5"
mpiexec -n 11 parallel 100 10 10 5 4 10
mpiexec -n 11 parallel 100 20 10 5 4 10
mpiexec -n 11 parallel 100 30 10 5 4 10
mpiexec -n 11 parallel 100 40 10 5 4 10
mpiexec -n 11 parallel 100 50 10 5 4 10
mpiexec -n 11 parallel 100 60 10 5 4 10
mpiexec -n 11 parallel 100 70 10 5 4 10
mpiexec -n 11 parallel 100 80 10 5 4 10

echo "10 x 10"
mpiexec -n 11 parallel 100 10 10 10 4 10
mpiexec -n 11 parallel 100 20 10 10 4 10
mpiexec -n 11 parallel 100 30 10 10 4 10
mpiexec -n 11 parallel 100 40 10 10 4 10
mpiexec -n 11 parallel 100 50 10 10 4 10
mpiexec -n 11 parallel 100 60 10 10 4 10
mpiexec -n 11 parallel 100 70 10 10 4 10
mpiexec -n 11 parallel 100 80 10 10 4 10

echo "10 x 20"
mpiexec -n 11 parallel 100 10 10 20 4 10
mpiexec -n 11 parallel 100 20 10 20 4 10
mpiexec -n 11 parallel 100 30 10 20 4 10
mpiexec -n 11 parallel 100 40 10 20 4 10
mpiexec -n 11 parallel 100 50 10 20 4 10
mpiexec -n 11 parallel 100 60 10 20 4 10
mpiexec -n 11 parallel 100 70 10 20 4 10
mpiexec -n 11 parallel 100 80 10 20 4 10

