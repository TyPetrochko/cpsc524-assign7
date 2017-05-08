#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

typedef struct chromosome {
  int *cover;
} chromosome;

typedef struct graph {
  int **adj_matrix;
} graph;

typedef struct coord {
  int x;
  int y;
} coord;

graph copyGraph(int n);

graph getRandomGraph(int n, int m);

chromosome getBlankChromosome(int n);

chromosome getRandomChromosome(int n);

graph getRandomGraph(int n, int m);

int evaluateFitness(chromosome c);

chromosome randomSolution();

chromosome crossover(chromosome c);

// return an array of four neighboring chromosomes
void getNeighbors(chromosome *buf, int x, int y, chromsome **torus, int width, int height);

// get coordinates of neighboring chromosomes
void getNeighborCoords(coord *buf, int x, int y, int width, int height);

// (in place) sort by fitness
void sortChromosomes(chromosome *chroms);

