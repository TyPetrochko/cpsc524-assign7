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

#define ROOT 0

#define SERIAL_ITER (10)
#define SERIAL_CHROMOSOMES (10) // must be divisible by two

#define TORUS_WIDTH (10)
#define TORUS_HEIGHT (6)

#define PARALLEL_ITERATIONS (100)
#define NUM_NEIGHBORS (4)

// chromosome serialGaTest(int n, int m, graph g){
// 
//   chromosome chroms[SERIAL_CHROMOSOMES];
//   int fitness[SERIAL_CHROMOSOMES];
// 
//   // init
//   for(int i = 0; i < SERIAL_CHROMOSOMES; i++){
//     chroms[i] = getRandomChromosome();
//     fitness[i] = evaluateFitness(chroms[i], g);
//   }
// 
//   for(int i = 0; i < SERIAL_ITER; i++){
//     sortChromosomes(chroms, g);
// 
//     // replace less-fit half with children from more-fit half
//     for(int j = SERIAL_CHROMOSOMES / 2; j < SERIAL_CHROMOSOMES; j++){
//       chroms[j] = 
//         crossover(
//             chroms[rand() % (SERIAL_CHROMOSOMES / 2)],
//             chroms[rand() % (SERIAL_CHROMOSOMES / 2)]);
//     }
// 
//     // re-eval fitness
//     for(int j = 0; j < SERIAL_CHROMOSOMES; j++){
//       fitness[j] = evaluateFitness(chroms[j], g);
//     }
//   }
// 
//   chromosome toReturn = chroms[0];
// 
//   // return optimal
//   for(int i = 0; i < SERIAL_CHROMOSOMES; i++){
//     if(evaluateFitness(chroms[i], g) > evaluateFitness(toReturn, g)){
//       toReturn = chroms[i];
//     }
//   }
// 
//   return toReturn;
// }

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
  
  printf("Random fitness eval: %d\n\n", evaluateFitness(c, g));

  // ========== BEGIN APPROX TEST ==========
  printf("STARTING APPROX TEST\n");
  
  timing(&wctime, &cputime);
  chromosome approx_solution = randomSolution(g);
  timing(&wctime_end, &cputime); 
  
  printf("\tSolution fitness eval: %d\n", evaluateFitness(approx_solution, g)); 
  printf("\tApprox time: %lf\n\n", wctime_end - wctime);
  
  // ========== BEGIN SERIAL TEST ==========
  printf("STARTING SERIAL TEST\n");
  
  timing(&wctime, &cputime);
  chromosome serial_solution = serialTorusTest(n, m, g);
  timing(&wctime_end, &cputime);
  
  printf("\tSolution fitness eval: %d\n", evaluateFitness(serial_solution, g));
  printf("\tSerial time: %lf\n\n", wctime_end - wctime);

  return 0;
}


