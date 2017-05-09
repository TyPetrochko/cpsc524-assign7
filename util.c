#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <omp.h>
#include "util.h"

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

int evaluateFitness(chromosome c, graph g){
  int sum = 0;
  for(int i = 0; i < g.n; i++){
    sum += c.cover[i];

    for(int j = i; j < g.n; j++){
      sum += g.n * (1 - c.cover[i]) * (1 - c.cover[j] * g.adj_matrix[i][j]);
    }
  }

  return (sum*(-1)); // we want to MINIMIZE this value
}

// TODO do we need this?
chromosome randomSolution(graph g){
  return getRandomChromosome(g.n);
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

// NOTE: Quicksort impl. derived from:
//    http://www.comp.dit.ie/rlawlor/Alg_DS/sorting/quickSort.c

int partition(chromosome a[], int l, int r, graph g);

// aux helper for quicksort
void quickSort(chromosome a[], int l, int r, graph g)
{
  int j;

  if( l < r ) 
  {
   // divide and conquer
   j = partition(a, l, r, g);
   quickSort(a, l, j-1, g);
   quickSort(a, j+1, r, g);
  }
}

// another aux helper for quicksort
int partition(chromosome a[], int l, int r, graph g) {
  int pivot, i, j;
  chromosome t;
  pivot = evaluateFitness(a[l], g);
  i = l; j = r+1;
   
  while(1)
  {
    do ++i; while( i <= r && evaluateFitness(a[i], g) <= pivot);
    do --j; while( evaluateFitness(a[j], g) > pivot );
    if( i >= j ) break;
    t = a[i]; a[i] = a[j]; a[j] = t;
  }
  t = a[l]; a[l] = a[j]; a[j] = t;
  return j;
}

void sortChromosomes(chromosome *chroms, graph g, int n){
  quickSort(chroms, 0, 3, g);
}

