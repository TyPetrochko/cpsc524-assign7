#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <omp.h>
#include "util.h"

#define MUTATE (1)

graph getRandomGraph(int n, int m){
  graph toReturn;

  toReturn.m = m;
  toReturn.n = n;
  toReturn.adj_matrix = (int**) malloc(sizeof(int*) * n);

  for(int i = 0; i < n; i++){
    toReturn.adj_matrix[i] = (int*) calloc(n, sizeof(int));
  }

  while(m > 0){
    int v = rand() % n;
    int w = rand() % n;

    if(v == w)
      continue;
    else if(toReturn.adj_matrix[v][w])
      continue;

    toReturn.adj_matrix[v][w] = 1;
    toReturn.adj_matrix[w][v] = 1;

    m--;
  }

  return toReturn;
}

graph copyGraph(graph g){
  graph toReturn;
  toReturn.n = g.n;
  toReturn.m = g.m;

  toReturn.adj_matrix = (int **) malloc(g.n * sizeof(int *));
  for(int i = 0; i < g.n; i++){
    toReturn.adj_matrix[i] = (int *) malloc(g.n * sizeof(int));
    for(int j = 0; j < g.n; j++){
      toReturn.adj_matrix[i][j] = g.adj_matrix[i][j];
    }
  }

  return toReturn;
}

chromosome getRandomChromosome(int n){
  chromosome toReturn;
  toReturn.n = n;
  toReturn.cover = (int*) malloc(sizeof(int) * n);

  for(int i = 0; i < n; i++){
    if(rand() % 2)
      toReturn.cover[i] = 1;
    else
      toReturn.cover[i] = 0;
  }

  return toReturn;
}

long evaluateFitness(chromosome c, graph g){
  long sum = 0;
  for(int i = 0; i < g.n; i++){
    sum += c.cover[i];

    if(c.cover[i] < 0 || c.cover[i] > 1){
      printf("PROGRAM ERROR: Bad chromosome!\n");
      exit(-1);
    }

    for(int j = i; j < g.n; j++){
      sum += g.n * (1 - c.cover[i]) * ((1 - c.cover[j]) * g.adj_matrix[i][j]);
    }
  }
  long a = (((g.n*g.n*g.n) - sum));
  if(sum < 0 || a < 0){
    printf("Yep, negative\n");
  }
  return a; // we want to MINIMIZE this value
}

// heuristic approach
chromosome randomSolution(graph g){
  chromosome toReturn;
  toReturn.n = g.n;
  toReturn.cover = (int *) calloc(g.n, sizeof(int));

  graph gprime = copyGraph(g);

  while(gprime.m > 0){
    // randomly choose an edge
    int i = 0;
    int j = 0;

    for(i = 0; i < g.n; i++){
      for(j = 0; j < g.n; j++){
        if(gprime.adj_matrix[i][j])
          goto escape;
      }
    }
escape:

    // remove the edge
    gprime.m--;
    gprime.adj_matrix[i][j] = 0;
    gprime.adj_matrix[j][i] = 0;

    // add the two vertices to the cover
    toReturn.cover[i] = 1;
    toReturn.cover[j] = 1;

    // remove all edges covered by vertices i and j
    for(int s = 0; s < g.n; s++){
      if(gprime.adj_matrix[i][s]){
        gprime.m--;
        gprime.adj_matrix[i][s] = 0;
        gprime.adj_matrix[s][i] = 0;
      }
      
      if(gprime.adj_matrix[j][s]){
        gprime.m--;
        gprime.adj_matrix[j][s] = 0;
        gprime.adj_matrix[s][j] = 0;
      }
    }
  }

  return toReturn;
}

void maybeMutate(chromosome c, int num, int denom){
  if(!MUTATE)
    return;

  if((rand() % denom) > (denom - num)){
    int flip1 = rand() % c.n;
    int flip2 = rand() % c.n;

    c.cover[flip1] = !c.cover[flip1];
    c.cover[flip2] = !c.cover[flip2];
  }
}

chromosome crossover(chromosome c, chromosome d){
  if(c.n != d.n){
    printf("Error: Incompatible chromosome length: %d vs %d", c.n, d.n);
    exit(-1);
  }
  
  chromosome toReturn = getRandomChromosome(c.n);

  int cross = rand() % c.n;

  for(int i = 0; i < c.n; i++){
    if(i < cross) toReturn.cover[i] = c.cover[i];
    else toReturn.cover[i] = d.cover[i];
  }

  return toReturn;
}

// return an array of four neighboring chromosomes
void getNeighbors(chromosome* buf, int x, int y, chromosome **torus, int width, int height){
  buf[0] = torus[(height + y - 1) % height][x]; // up
  buf[1] = torus[(y + 1) % height][x]; // down
  buf[2] = torus[y][(x + 1) % width]; // right
  buf[3] = torus[y][(x + width - 1) % width]; // left
}

// only need the top two sorted...
void sortChromosomes(chromosome *chroms, graph g, int n){
  for(int i = 0; i < 1; i++){
    if(evaluateFitness(chroms[i], g) > evaluateFitness(chroms[3], g)){
      chromosome tmp = chroms[i];
      chroms[i] = chroms[3];
      chroms[3] = tmp;
    } else if(evaluateFitness(chroms[i], g) > evaluateFitness(chroms[2], g)){
      chromosome tmp = chroms[i];
      chroms[i] = chroms[2];
      chroms[2] = tmp;
    }
  }
}

