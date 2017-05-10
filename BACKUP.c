#include "timing.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h> 
#include <stdbool.h>
#include <math.h>
#include <omp.h>
#include "mpi.h"
#include "util.h"

#define ROOT (0)
#define TAG (99)

#define SERIAL_ITER (10)
#define SERIAL_CHROMOSOMES (10) // must be divisible by two

#define PARALLEL_ITERATIONS (100)
#define NUM_NEIGHBORS (4)

int TORUS_WIDTH;
int TORUS_HEIGHT;

int NUMERATOR;
int DENOMINATOR;

chromosome serialTorusTest(int n, int m, graph g){
  // we need a backup torus to copy new values into
  chromosome **torus;
  chromosome **backup;

  torus  = (chromosome **) malloc(sizeof(chromosome) * TORUS_HEIGHT);
  backup = (chromosome **) malloc(sizeof(chromosome) * TORUS_HEIGHT);

  // initialize torus
  for(int i = 0; i < TORUS_HEIGHT; i++){
    torus[i]  = (chromosome *) malloc(TORUS_WIDTH * sizeof(chromosome));
    backup[i] = (chromosome *) malloc(TORUS_WIDTH * sizeof(chromosome));
    for(int j = 0; j < TORUS_WIDTH; j++){
      torus[i][j] = getRandomChromosome(n);
    }
  }

  // main processing loop
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
          maybeMutate(backup[i][j], NUMERATOR, DENOMINATOR);
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

int main(int argc, char **argv){
  // ============= BEGIN MAIN =============
  srand(24); // lucky number

  if(argc != 7){
    printf(" \n");
    printf(" Usage: ./parallel [n] [m] [width] [height] [num] [denom] \n");
    printf(" \n");
    printf("     - num/denom is the mutation likelihood\n");
    printf("     - width and height are the torus dimensions\n");
    printf(" \n");
    printf(" \n");
    exit(-1);
  }

  int n, m;

  n = atoi(argv[1]);
  m = atoi(argv[2]);
  TORUS_WIDTH = atoi(argv[3]);
  TORUS_HEIGHT = atoi(argv[4]);
  NUMERATOR = atoi(argv[5]);
  DENOMINATOR = atoi(argv[6]);


  chromosome c = getRandomChromosome(n);
  
  double wctime, wctime_end, cputime;
  
  // ============= BEGIN PARALLEL CODE =============
  int size, rank;
  MPI_Status status;

  MPI_Init(&argc,&argv);
  
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  srand(rank);

  if(rank == ROOT){
    timing(&wctime, &cputime);
    
    graph g = getRandomGraph(n, m);

    // send the graph to all workers
    for(int i = 0; i < n; i++){
      MPI_Bcast(g.adj_matrix[i], n, MPI_INT, ROOT, MPI_COMM_WORLD);
    }
   
    // let every iteration happen
    for(int iter = 0; iter < PARALLEL_ITERATIONS; iter++){
      MPI_Barrier(MPI_COMM_WORLD);
    }

    // initialize a new chromosome (for buffer space)
    chromosome solution = getRandomChromosome(n);
    chromosome next = getRandomChromosome(n);

    // receive all chromosomes from all workers
    for(int i = 0; i < TORUS_HEIGHT * TORUS_WIDTH; i++){
      MPI_Recv(
          next.cover, 
          n,
          MPI_INT,
          MPI_ANY_SOURCE,
          MPI_ANY_TAG,
          MPI_COMM_WORLD,
          &status);

      if(evaluateFitness(next, g) > evaluateFitness(solution, g))
        solution = next;
    }

    long parallel_time, parallel_fitness;
    long approx_time, approx_fitness;
    long serial_time, serial_fitness;
    long random_fitness = evaluateFitness(getRandomChromosome(n), g);
    
    timing(&wctime_end, &cputime);
    parallel_time = wctime_end - wctime;
    parallel_fitness = evaluateFitness(solution, g) - random_fitness;
    MPI_Finalize();
    
    // ========== BEGIN APPROX TEST ==========
    
    timing(&wctime, &cputime);
    chromosome approx_solution = randomSolution(g);
    timing(&wctime_end, &cputime); 
    
    approx_time = wctime_end - wctime;
    approx_fitness = evaluateFitness(approx_solution, g) - random_fitness;
    
    // ========== BEGIN SERIAL TEST ==========
    
    timing(&wctime, &cputime);
    chromosome serial_solution = serialTorusTest(n, m, g);
    timing(&wctime_end, &cputime);
    
    serial_time = wctime_end - wctime;
    serial_fitness = evaluateFitness(serial_solution, g) - random_fitness;
   
    printf("%8lf \t %8lf\t", serial_time, serial_fitness);
    printf("%8lf \t %8lf\t", approx_time, approx_fitness);
    printf("%8lf \t %8lf\t\n", parallel_time, parallel_fitness);
    return 0;
  } else {
    // WORKER CODE
    
    // receive graph
    graph g = getRandomGraph(n, m); // does mallocing for us
    g.m = m;
    g.n = n;

    for(int i = 0; i < n; i++){
      MPI_Bcast(g.adj_matrix[i], n, MPI_INT, ROOT, MPI_COMM_WORLD);
    }
   
    chromosome *chrom = malloc(TORUS_HEIGHT * sizeof(chromosome));
    chromosome *back  = malloc(TORUS_HEIGHT * sizeof(chromosome));

    // initialize random chromosomes
    for(int i = 0; i < TORUS_HEIGHT; i++){
      chrom[i] = getRandomChromosome(n);
    }

    int worker_idx = rank-1;
      
    chromosome left[TORUS_HEIGHT];
    chromosome right[TORUS_HEIGHT];

    // handle mallocing for us!
    for(int i = 0; i < TORUS_HEIGHT; i++){
      left[i] = getRandomChromosome(n);
      right[i] = getRandomChromosome(n);
    }
    
    for(int iter = 0; iter < PARALLEL_ITERATIONS; iter++){
      MPI_Barrier(MPI_COMM_WORLD);

    
      // ============= BEGIN GENETIC EXCHANGE =============
      // ============= EVEN WORKERS =============
      if(worker_idx % 2 == 0){
        // send info to right
        for(int i = 0; i < TORUS_HEIGHT; i++){
          MPI_Send(
              chrom[i].cover, 
              n,
              MPI_INT,
              (worker_idx + 1) % TORUS_WIDTH + 1,
              TAG,
              MPI_COMM_WORLD);
        }
        // send info to left
        for(int i = 0; i < TORUS_HEIGHT; i++){
          MPI_Send(
              chrom[i].cover, 
              n,
              MPI_INT,
              (worker_idx + TORUS_WIDTH - 1) % TORUS_WIDTH + 1,
              TAG,
              MPI_COMM_WORLD);
        }
        // receive info from left
        for(int i = 0; i < TORUS_HEIGHT; i++){
          MPI_Recv(
              left[i].cover, 
              n,
              MPI_INT,
              (worker_idx + TORUS_WIDTH - 1) % TORUS_WIDTH + 1,
              MPI_ANY_TAG,
              MPI_COMM_WORLD,
              &status);
        }
        // receive info from right
        for(int i = 0; i < TORUS_HEIGHT; i++){
          MPI_Recv(
              right[i].cover, 
              n,
              MPI_INT,
              (worker_idx + 1) % TORUS_WIDTH + 1,
              MPI_ANY_TAG,
              MPI_COMM_WORLD,
              &status);
        }
        // ============= ODD WORKERS =============
      } else{
        // receive info from left
        for(int i = 0; i < TORUS_HEIGHT; i++){
          MPI_Recv(
              left[i].cover, 
              n,
              MPI_INT,
              (worker_idx + TORUS_WIDTH - 1) % TORUS_WIDTH + 1,
              MPI_ANY_TAG,
              MPI_COMM_WORLD,
              &status);
        }
        // receive info from right
        for(int i = 0; i < TORUS_HEIGHT; i++){
          MPI_Recv(
              right[i].cover, 
              n,
              MPI_INT,
              (worker_idx + 1) % TORUS_WIDTH + 1,
              MPI_ANY_TAG,
              MPI_COMM_WORLD,
              &status);
        }
        // send info to right
        for(int i = 0; i < TORUS_HEIGHT; i++){
          MPI_Send(
              chrom[i].cover, 
              n,
              MPI_INT,
              (worker_idx + 1) % TORUS_WIDTH + 1,
              TAG,
              MPI_COMM_WORLD);
        }
        // send info to left
        for(int i = 0; i < TORUS_HEIGHT; i++){
          MPI_Send(
              chrom[i].cover, 
              n,
              MPI_INT,
              (worker_idx + TORUS_WIDTH - 1) % TORUS_WIDTH + 1,
              TAG,
              MPI_COMM_WORLD);
        }
      }
      // ============= END GENETIC EXCHANGE =============
      
      // mate (unless we're the most superior
      for(int i = 0; i < TORUS_HEIGHT; i++){
        chromosome neighbors[4];

        neighbors[0] = chrom[(i + 1) % TORUS_HEIGHT];
        neighbors[1] = chrom[(i + TORUS_HEIGHT - 1) % TORUS_HEIGHT];
        neighbors[2] = left[i];
        neighbors[3] = right[i];
        sortChromosomes(neighbors, g, 4);
       
        // find genetically superior mates
        if(evaluateFitness(chrom[i], g) >= evaluateFitness(neighbors[3], g)){
          back[i] = chrom[i];
          maybeMutate(back[i], NUMERATOR, DENOMINATOR);
        } else {
          back[i] = crossover(neighbors[2], neighbors[3]);
        }
      }

      // switch buffers
      chromosome *tmp = chrom;
      chrom = back;
      back = tmp;
    
    } // end iter loop

    // send everything to master
    for(int i = 0; i < TORUS_HEIGHT; i++){
      MPI_Send(
          chrom[i].cover, 
          n,
          MPI_INT,
          ROOT,
          TAG,
          MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
  } // end if-master-worker conditional
} 

