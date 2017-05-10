module load Langs/Intel/15
module load Langs/Intel/15 MPI/OpenMPI/1.8.6-intel15

echo "format: serial/approx/parallel time/fitness" echo "10 x 5"
# mpiexec -n 11 parallel 200 100 10 5 5 10
# mpiexec -n 11 parallel 400 200 10 5 5 10
# mpiexec -n 11 parallel 600 300 10 5 5 10
# mpiexec -n 11 parallel 800 400 10 5 5 10
# mpiexec -n 11 parallel 1000 500 10 5 5 10
# mpiexec -n 11 parallel 1200 600 10 5 5 10
# mpiexec -n 11 parallel 1400 700 10 5 5 10
# mpiexec -n 11 parallel 1600 700 10 5 5 10
# mpiexec -n 11 parallel 1800 700 10 5 5 10
# mpiexec -n 11 parallel 2000 700 10 5 5 10
# mpiexec -n 11 parallel 2200 700 10 5 5 10
# mpiexec -n 11 parallel 2400 700 10 5 5 10
mpiexec -n 11 parallel 2600 700 10 5 5 10
mpiexec -n 11 parallel 2800 700 10 5 5 10
mpiexec -n 11 parallel 3000 700 10 5 5 10

