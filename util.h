#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

typedef struct chromosome {
  int n;
  int *cover;
} chromosome;

typedef struct graph {
  int **adj_matrix;
  int n;
  int m;
} graph;

typedef struct coord {
  int x;
  int y;
} coord;

graph getRandomGraph(int n, int m);

graph copyGraph(graph g);

chromosome getRandomChromosome(int n);

graph getRandomGraph(int n, int m);

int evaluateFitness(chromosome c, graph g);

// TODO do we need this?
chromosome randomSolution(graph g);

chromosome crossover(chromosome c, chromosome d);

// return an array of four neighboring chromosomes
void getNeighbors(chromosome *buf, int x, int y, chromosome **torus, int width, int height);

// DO NOT NEED THIS FOR NOW
// get coordinates of neighboring chromosomes
// void getNeighborCoords(coord *buf, int x, int y, int width, int height);

// (in place) sort by fitness
void sortChromosomes(chromosome *chroms, graph g, int n);

