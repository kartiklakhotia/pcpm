#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#define GRAPH_HEADER_INCL

#if defined (HUGE_EDGE) || defined (HUGE_VERTEX)
typedef unsigned long long int intE;
#define LSB_MASK 0x0000000000000001
#else
#define LSB_MASK 0x00000001
typedef unsigned int intE;
#endif
//typedef unsigned long long int intE;
//typedef unsigned int intV;

#ifdef HUGE_VERTEX
typedef unsigned long long int intV;
#else
typedef unsigned int intV;
#endif


typedef struct graph 
{
    intV numVertex;
    intE numEdges;
    intE* VI;
    intV* EI;
    unsigned int* EW;
    bool weighted;
} graph;

int read_csr (char*, graph*);

void printGraph (graph*);

void write_csr (char*, graph*);

void createReverseCSR (graph*, graph*);

void freeMem(graph*);

void sortEdges(graph*);
