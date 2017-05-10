#include "timing.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h> 
#include <stdbool.h>
#include <math.h>
#include <omp.h>
#include "util.h"

#define ROOT (0)
#define TAG (99)

#define SERIAL_ITER (10)
#define SERIAL_CHROMOSOMES (10) // must be divisible by two

#define TORUS_WIDTH (10)
#define TORUS_HEIGHT (6)

#define PARALLEL_ITERATIONS (100)
#define NUM_NEIGHBORS (4)

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
        if(evaluateFitness(torus[i][j], g) >= evaluateFitness(neighbors[NUM_NEIGHBORS - 1], g)){
          backup[i][j] = torus[i][j];
          maybeMutate(backup[i][j]);
        }else{
          backup[i][j] = crossover(neighbors[NUM_NEIGHBORS - 1], neighbors[NUM_NEIGHBORS - 2]);
        }
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

// int main(int argc, char **argv){
//   // ============= BEGIN PARALLEL CODE =============
//   int size, rank;
//   MPI_Status status;
// 
//   MPI_Init(&argc,&argv);
//   
//   MPI_Comm_size(MPI_COMM_WORLD,&size);
//   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//   
//   srand(rank); // lucky number
// 
//   if(rank == ROOT){
//     printf("STARTING PARALLEL TEST:\n");
//     timing(&wctime, &cputime);
//     
//     graph g = getRandomGraph(n, m);
// 
//     for(int i = 0; i < n; i++){
//       MPI_Bcast(g.adj_matrix[i], n, MPI_INT, ROOT, MPI_COMM_WORLD);
//     }
//     
//     for(int iter = 0; iter < PARALLEL_ITERATIONS; iter++){
//       MPI_Barrier(MPI_COMM_WORLD);
//     }
// 
// 
//     chromosome solution = getRandomChromosome(n);
//     chromosome next = getRandomChromosome(n);
// 
//     for(int i = 0; i < TORUS_HEIGHT * TORUS_WIDTH; i++){
//       MPI_Recv(
//           next.cover, 
//           n,
//           MPI_INT,
//           MPI_ANY_SOURCE,
//           MPI_ANY_TAG,
//           MPI_COMM_WORLD,
//           &status);
// 
//       if(evaluateFitness(next, g) > evaluateFitness(solution, g))
//         solution = next;
//     }
//     
//     timing(&wctime_end, &cputime);
//     printf("\tParallel time: %lf\n", wctime_end - wctime);
//     printf("\tParallel fitness eval: %lf\n\n", evaluateFitness(solution, g));
//     MPI_Finalize();
//     
//     // ========== BEGIN APPROX TEST ==========
//     printf("STARTING APPROX TEST\n");
//     
//     timing(&wctime, &cputime);
//     chromosome approx_solution = randomSolution(g);
//     timing(&wctime_end, &cputime); 
//     
//     printf("\tApprox fitness eval: %lf\n", evaluateFitness(approx_solution, g)); 
//     printf("\tApprox time: %lf\n\n", wctime_end - wctime);
//     
//   } else {
//     // WORKER CODE
//     
//     // receive graph
//     graph g = getRandomGraph(n, m); // does mallocing for us
//     g.m = m;
//     g.n = n;
// 
//     for(int i = 0; i < n; i++){
//       MPI_Bcast(g.adj_matrix[i], n, MPI_INT, ROOT, MPI_COMM_WORLD);
//     }
//    
//     // chromosome chrom[TORUS_HEIGHT];
//     // chromosome back[TORUS_HEIGHT];
//     chromosome *chrom = malloc(TORUS_HEIGHT * sizeof(chromosome));
//     chromosome *back  = malloc(TORUS_HEIGHT * sizeof(chromosome));
// 
//     for(int i = 0; i < TORUS_HEIGHT; i++){
//       chrom[i] = getRandomChromosome(n);
//     }
// 
//     int worker_idx = rank-1;
//       
//     chromosome left[TORUS_HEIGHT];
//     chromosome right[TORUS_HEIGHT];
// 
//     // handle mallocing for us!
//     for(int i = 0; i < TORUS_HEIGHT; i++){
//       left[i] = getRandomChromosome(n);
//       right[i] = getRandomChromosome(n);
//     }
//     
//     for(int iter = 0; iter < PARALLEL_ITERATIONS; iter++){
//       MPI_Barrier(MPI_COMM_WORLD);
// 
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
//               (worker_idx + 1) % TORUS_WIDTH + 1,
//               TAG,
//               MPI_COMM_WORLD);
//         }
//         // send info to left
//         for(int i = 0; i < TORUS_HEIGHT; i++){
//           MPI_Send(
//               chrom[i].cover, 
//               n,
//               MPI_INT,
//               (worker_idx + TORUS_WIDTH - 1) % TORUS_WIDTH + 1,
//               TAG,
//               MPI_COMM_WORLD);
//         }
//         // receive info from left
//         for(int i = 0; i < TORUS_HEIGHT; i++){
//           MPI_Recv(
//               left[i].cover, 
//               n,
//               MPI_INT,
//               (worker_idx + TORUS_WIDTH - 1) % TORUS_WIDTH + 1,
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
//               (worker_idx + 1) % TORUS_WIDTH + 1,
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
//               (worker_idx + TORUS_WIDTH - 1) % TORUS_WIDTH + 1,
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
//               (worker_idx + 1) % TORUS_WIDTH + 1,
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
//               (worker_idx + 1) % TORUS_WIDTH + 1,
//               TAG,
//               MPI_COMM_WORLD);
//         }
//         // send info to left
//         for(int i = 0; i < TORUS_HEIGHT; i++){
//           MPI_Send(
//               chrom[i].cover, 
//               n,
//               MPI_INT,
//               (worker_idx + TORUS_WIDTH - 1) % TORUS_WIDTH + 1,
//               TAG,
//               MPI_COMM_WORLD);
//         }
//       }
//       // ============= END GENETIC EXCHANGE =============
//       
//       // mate (unless we're the most superior
//       for(int i = 0; i < TORUS_HEIGHT; i++){
//         chromosome neighbors[4];
//         
// 
//         neighbors[0] = chrom[(i + 1) % TORUS_HEIGHT];
//         neighbors[1] = chrom[(i + TORUS_HEIGHT - 1) % TORUS_HEIGHT];
//         neighbors[2] = left[i];
//         neighbors[3] = right[i];
//         sortChromosomes(neighbors, g, 4);
//         
//         if(evaluateFitness(chrom[i], g) >= evaluateFitness(neighbors[3], g)){
//           back[i] = chrom[i];
//           maybeMutate(back[i]);
//         } else {
//           back[i] = crossover(neighbors[2], neighbors[3]);
//         }
//       }
// 
//       // switch buffers
//       chromosome *tmp = chrom;
//       chrom = back;
//       back = tmp;
//     
//     } // end iter loop
// 
//     for(int i = 0; i < TORUS_HEIGHT; i++){
//       MPI_Send(
//           chrom[i].cover, 
//           n,
//           MPI_INT,
//           ROOT,
//           TAG,
//           MPI_COMM_WORLD);
//     }
// 
//     MPI_Finalize();
//     return 0;
//   } // end if-master-worker conditional
// } 

int main(int argc, char**argv){
  // ============= BEGIN MAIN =============
  srand(24); // lucky number

  if(argc != 3){
    printf("Usage: ./parallel [n] [m]\n");
    exit(-1);
  }

  int n, m;

  n = atoi(argv[1]);
  m = atoi(argv[2]);

  chromosome c = getRandomChromosome(n);
  
  double wctime, wctime_end, cputime;
    
  graph g = getRandomGraph(n, m);
  
  // ========== BEGIN SERIAL TEST ==========
  printf("STARTING SERIAL TEST\n");
  
  timing(&wctime, &cputime);
  chromosome serial_solution = serialTorusTest(n, m, g);
  timing(&wctime_end, &cputime);
  
  printf("\tSerial Solution fitness eval: %lf\n", evaluateFitness(serial_solution, g));
  printf("\tSerial time: %lf\n\n", wctime_end - wctime);
  return 0;
}

