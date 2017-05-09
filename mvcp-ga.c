#include "timing.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h> 
#include <stdbool.h>
#include <math.h>
#include <omp.h>
// #include "mpi.h"
#include "util.h"

#define ROOT 0

#define SERIAL_ITER (10)
#define SERIAL_CHROMOSOMES (10) // must be divisible by two

#define TORUS_WIDTH (20)
#define TORUS_HEIGHT (10)

#define PARALLEL_ITERATIONS (10)
#define NUM_NEIGHBORS (4)

/*
chromosome serialGaTest(int n, int m, graph g){

  chromosome chroms[SERIAL_CHROMOSOMES];
  int fitness[SERIAL_CHROMOSOMES];

  // init
  for(int i = 0; i < SERIAL_CHROMOSOMES; i++){
    chroms[i] = getRandomChromosome();
    fitness[i] = evaluateFitness(chroms[i], g);
  }

  for(int i = 0; i < SERIAL_ITER; i++){
    sortChromosomes(chroms, g);

    // replace less-fit half with children from more-fit half
    for(int j = SERIAL_CHROMOSOMES / 2; j < SERIAL_CHROMOSOMES; j++){
      chroms[j] = 
        crossover(
            chroms[rand() % (SERIAL_CHROMOSOMES / 2)],
            chroms[rand() % (SERIAL_CHROMOSOMES / 2)]);
    }

    // re-eval fitness
    for(int j = 0; j < SERIAL_CHROMOSOMES; j++){
      fitness[j] = evaluateFitness(chroms[j], g);
    }
  }

  chromosome toReturn = chroms[0];

  // return optimal
  for(int i = 0; i < SERIAL_CHROMOSOMES; i++){
    if(evaluateFitness(chroms[i], g) > evaluateFitness(toReturn, g)){
      toReturn = chroms[i];
    }
  }

  return toReturn;
}
*/


chromosome serialTorusTest(int n, int m, graph g){
 
  // we need a backup torus to copy new values into
  chromosome **torus;
  chromosome **backup;

  torus  = (chromosome **) malloc(sizeof(chromosome) * TORUS_HEIGHT);
  backup = (chromosome **) malloc(sizeof(chromosome) * TORUS_HEIGHT);

  for(int i = 0; i < TORUS_HEIGHT; i++){
    torus[i]  = (chromosome *) malloc(TORUS_WIDTH * sizeof(chromosome));
    backup[i] = (chromosome *) malloc(TORUS_WIDTH * sizeof(chromosome));
    for(int j = 0; j < TORUS_WIDTH; j++){
      torus[i][j] = getRandomChromosome(n);
    }

    // sortChromosomes(torus[i], g);
  }

  for(int iter = 0; iter < PARALLEL_ITERATIONS; iter++){

    for(int i = 0; i < TORUS_HEIGHT; i++){
      for(int j = 0; j < TORUS_WIDTH; j++){

        // get all 4 neighboring chromosomes
        chromosome neighbors[NUM_NEIGHBORS];
        getNeighbors(neighbors, j, i, torus, TORUS_WIDTH, TORUS_HEIGHT);

        sortChromosomes(neighbors, g, NUM_NEIGHBORS);

        // replace this chrom with a crossover of fittest neighbors (or don't)
        if(evaluateFitness(torus[i][j], g) > evaluateFitness(neighbors[0], g))
          backup[i][j] = torus[i][j];
        else
          backup[i][j] = crossover(neighbors[0], neighbors[1]);
      }
    } // end i-j loop

    // switch the buffers
    chromosome **tmp = torus;
    torus = backup;
    backup = tmp;
  } // end iter loop
  
  // return optimal chromosome
  chromosome toReturn = torus[0][0];
  for(int i = 0; i < TORUS_HEIGHT; i++){
    for(int j = 0; j < TORUS_WIDTH; j++){
      if(evaluateFitness(torus[i][j], g) > evaluateFitness(toReturn, g)){
        toReturn = torus[i][j];
      }
    }
  }

  return toReturn;
}

// chromosome mpiGaTest(int n, int m, int argc, char **argv, graph g){
//   int size, rank;
//   MPI_Status status;
// 
//   MPI_Init(&argc,&argv);
//   
//   MPI_Comm_size(MPI_COMM_WORLD,&size);
//   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
// 
//   if(rank == ROOT){
//     // MASTER CODE
//     for(int iter = 0; iter < PARALLEL_ITERATIONS; iter++){
//       MPI_Barrier(MPI_COMM_WORLD);
//       printf("Begin iteration %d\n", iter);
//     }
// 
//     chromosome toReturn = getRandomChromosome(n);
// 
//     for(int i = 0; i < TORUS_HEIGHT * TORUS_WIDTH; i++){
//       chromosome next;
//       MPI_Recv(
//           next.cover, 
//           n,
//           MPI_INT,
//           MPI_ANY_SOURCE,
//           MPI_ANY_TAG,
//           MPI_COMM_WORLD,
//           &status);
// 
//       if(evaluateFitness(next, g) > evaluateFitness(toReturn, g))
//         toReturn = next;
//     }
// 
//     return toReturn;
//   } else {
//     // WORKER CODE
//    
//     chromosome chrom[TORUS_HEIGHT];
// 
//     for(int i = 0; i < TORUS_HEIGHT; i++){
//       chrom[i] = getRandomChromosome(n);
//     }
// 
//     int worker_idx = rank-1;
//     
//     for(int iter = 0; iter < PARALLEL_ITERATIONS; iter++){
//       MPI_Barrier(MPI_COMM_WORLD);
// 
//       chromosome left[TORUS_HEIGHT];
//       chromosome right[TORUS_WIDTH];
//     
//       // ============= BEGIN GENETIC EXCHANGE =============
//       // ============= EVEN WORKERS =============
//       if(worker_idx % 2 == 0){
//         // send info to right
//         for(int i = 0; i < TORUS_HEIGHT; i++){
//           MPI_Send(
//               chrom[i].cover, 
//               n,
//               MPI_INT,
//               (worker_idx + 1) % TORUS_WIDTH,
//               MPI_ANY_TAG,
//               MPI_COMM_WORLD);
//         }
//         // send info to left
//         for(int i = 0; i < TORUS_HEIGHT; i++){
//           MPI_Send(
//               chrom[i].cover, 
//               n,
//               MPI_INT,
//               (worker_idx + worker_idx - 1) % TORUS_WIDTH,
//               MPI_ANY_TAG,
//               MPI_COMM_WORLD);
//         }
//         // receive info from left
//         for(int i = 0; i < TORUS_HEIGHT; i++){
//           MPI_Recv(
//               left[i].cover, 
//               n,
//               MPI_INT,
//               (worker_idx + worker_idx - 1) % TORUS_WIDTH,
//               MPI_ANY_TAG,
//               MPI_COMM_WORLD,
//               &status);
//         }
//         // receive info from right
//         for(int i = 0; i < TORUS_HEIGHT; i++){
//           MPI_Recv(
//               right[i].cover, 
//               n,
//               MPI_INT,
//               (worker_idx + 1) % TORUS_WIDTH,
//               MPI_ANY_TAG,
//               MPI_COMM_WORLD,
//               &status);
//         }
//         // ============= ODD WORKERS =============
//       } else{
//         // receive info from left
//         for(int i = 0; i < TORUS_HEIGHT; i++){
//           MPI_Recv(
//               left[i].cover, 
//               n,
//               MPI_INT,
//               (worker_idx + worker_idx - 1) % TORUS_WIDTH,
//               MPI_ANY_TAG,
//               MPI_COMM_WORLD,
//               &status);
//         }
//         
//         // receive info from right
//         for(int i = 0; i < TORUS_HEIGHT; i++){
//           MPI_Recv(
//               right[i].cover, 
//               n,
//               MPI_INT,
//               (worker_idx + 1) % TORUS_WIDTH,
//               MPI_ANY_TAG,
//               MPI_COMM_WORLD,
//               &status);
//         }
//         // send info to right
//         for(int i = 0; i < TORUS_HEIGHT; i++){
//           MPI_Send(
//               chrom[i].cover, 
//               n,
//               MPI_INT,
//               (worker_idx + 1) % TORUS_WIDTH,
//               MPI_ANY_TAG,
//               MPI_COMM_WORLD);
//         }
//         // send info to left
//         for(int i = 0; i < TORUS_HEIGHT; i++){
//           MPI_Send(
//               chrom[i].cover, 
//               n,
//               MPI_INT,
//               (worker_idx + worker_idx - 1) % TORUS_WIDTH,
//               MPI_ANY_TAG,
//               MPI_COMM_WORLD);
//         }
//       }
//       // ============= END GENETIC EXCHANGE =============
//       
//       // mate (unless we're the most superior
//       for(int i = 0; i < TORUS_HEIGHT; i++){
//         chromosome neighbors[4];
// 
//         neighbors[0] = chrom[(i + 1) % TORUS_HEIGHT];
//         neighbors[1] = chrom[(i + TORUS_HEIGHT - 1) % TORUS_HEIGHT];
//         neighbors[2] = left[i];
//         neighbors[3] = right[i];
//         
//         sortChromosomes(neighbors, g);
//         if(evaluateFitness(chrom[i], g) > evaluateFitness(neighbors[0], g))
//           continue;
//         else
//           chrom[i] = crossover(neighbors[0], neighbors[1]);
//       }
//     } // end iter loop
// 
//     for(int i = 0; i < TORUS_HEIGHT; i++){
//       MPI_Send(
//           chrom[i].cover, 
//           n,
//           MPI_INT,
//           ROOT,
//           MPI_ANY_TAG,
//           MPI_COMM_WORLD);
//     }
//     exit(0);
//   } // end if-master-worker conditional
// } 

int main(int argc, char **argv){
  srand(24); // lucky number

  if(argc != 3){
    printf("Usage: ./mvcp-ga [n] [m]\n");
    exit(-1);
  }

  int n, m;

  n = atoi(argv[1]);
  m = atoi(argv[2]);

  graph g = getRandomGraph(n, m);
  chromosome c = getRandomChromosome(n);
  
  double wctime, wctime_end, cputime;

  printf("STARTING SERIAL TEST\n");
  timing(&wctime, &cputime);
  chromosome serial_solution = serialTorusTest(n, m, g);
  timing(&wctime_end, &cputime);
  
  printf("\tRandom fitness eval: %d\n", evaluateFitness(c, g));
  printf("\tSolution fitness eval: %d\n", evaluateFitness(serial_solution, g));
  
  printf("Serial time: %lf\n", wctime_end - wctime);

  // timing(&wctime, &cputime);
  // chromosome mpi_solution = mpiGaTest(n, m, argc, argv, g);
  // timing(&wctime_end, &cputime);
  // 
  // printf("Parallel time: %lf", wctime_end - wctime);
  return 0;
}

