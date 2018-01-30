#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#ifndef SORT_HEADER
#include "sort.hpp"
#endif

#define GRAPH_HEADER

typedef struct graph 
{
    unsigned int numVertex;
    unsigned int numEdges;
    unsigned int* VI;
    unsigned int* EI;
    float* attr;
} graph;

typedef struct edge
{
    unsigned int src;
    unsigned int dst;
} edge;

int read_csr (char* filename, graph* G)
{
    FILE* graphFile = fopen(filename, "rb");
    if (graphFile == NULL)
    {
        fputs("file error", stderr);
        return -1;
    }
    fread (&(G->numVertex), sizeof(unsigned int), 1, graphFile);
    
    fread (&(G->numEdges), sizeof(unsigned int), 1, graphFile);

    G->VI = new unsigned int[G->numVertex+1];
    fread (G->VI, sizeof(unsigned int), G->numVertex, graphFile);
    if (feof(graphFile))
    {
        delete[] G->VI;
        printf("unexpected end of file while reading vertices\n");
        return -1;
    }
    else if (ferror(graphFile))
    {
        delete[] G->VI;
        printf("error reading file\n");
        return -1;
    }
    G->VI[G->numVertex] = G->numEdges;

    G->EI = new unsigned int[G->numEdges];
    fread (G->EI, sizeof(unsigned int), G->numEdges, graphFile);
    if (feof(graphFile))
    {
        delete[] G->EI;
        delete[] G->VI;
        printf("unexpected end of file while reading edges\n");
        return -1;
    }
    else if (ferror(graphFile))
    {
        delete[] G->EI;
        delete[] G->VI;
        printf("error reading file\n");
        return -1;
    }
    fclose(graphFile);

    return 1;
}

void write_csr (char* filename, graph* G)
{
    FILE* fp = fopen(filename, "wb");
    if (fp == NULL)
    {
        fputs("file error", stderr);
        return;
    }
    fwrite(&G->numVertex, sizeof(unsigned int), 1, fp); 
    fwrite(&G->numEdges, sizeof(unsigned int), 1, fp); 
    fwrite(G->VI, sizeof(unsigned int), G->numVertex, fp); 
    fwrite(G->EI, sizeof(unsigned int), G->numEdges, fp); 
    fclose(fp); 
}

void printGraph(graph* G)
{
    printf("num vertices = %d\n numEdges = %d\n", G->numVertex, G->numEdges);
    for (unsigned int i=0; i<=G->numVertex; i++)
        printf("%d  ", G->VI[i]);
    printf("\n");
    for (unsigned int i=0; i<G->numEdges; i++)
        printf("%d  ", G->EI[i]);
    printf("\n"); 
}

void transposeCSR(graph* G1)
{
    unsigned int* newVI = new unsigned int[G1->numVertex+1]();
    unsigned int* newEI = new unsigned int[G1->numEdges]; 

    for (unsigned int i=0; i<G1->numEdges; i++)
    {
        newVI[G1->EI[i]+1]++;
    }
    for (unsigned int i=0; i<G1->numVertex; i++)
        newVI[i+1] += newVI[i];

    unsigned int* tempId = new unsigned int [G1->numVertex]();
    for (unsigned int i=0; i<G1->numVertex; i++)
    {
        for (unsigned int j=G1->VI[i]; j<G1->VI[i+1]; j++)
        {
            newEI[newVI[G1->EI[j]] + tempId[G1->EI[j]]] = i;
            tempId[G1->EI[j]]++;
        } 
    }
    delete[] G1->VI;
    delete[] G1->EI;
    delete[] tempId;
    G1->VI = newVI;
    G1->EI = newEI;
}

void sortEdges(graph* G)
{
    #pragma omp parallel for
    for (unsigned int i=0; i<G->numVertex; i++)
    {
        if (G->VI[i+1] > (G->VI[i]+1))
            mergeSortWOkey<unsigned int>(G->EI, G->VI[i], G->VI[i+1]-1);
    }
    return;
}

void findOutDeg(graph* G, unsigned int* degArr)
{
    #pragma omp parallel for
    for (unsigned int i=0; i<G->numVertex; i++)
        degArr[i] = G->VI[i+1] - G->VI[i];
    return;
}

void findOutDegTrans(graph* G, unsigned int* degArr)
{
    for (unsigned int i=0; i<G->numEdges; i++)
        degArr[G->EI[i]]++;
    return;
}

void initGraph (graph* G, unsigned int* outDeg)
{
    G->attr = new float [G->numVertex];
   #pragma omp parallel for
    for (unsigned int i=0; i<G->numVertex; i++)
    {
        G->attr[i]   = (1.0);  
    }
    return;
}

void initGraphPullPR (graph* G, unsigned int* outDeg)
{
    G->attr = new float [G->numVertex]();

    #pragma omp parallel for
    for (unsigned int i=0; i<G->numVertex; i++)
    {
        G->attr[i] = 1;
        if (outDeg[i] > 0)
            G->attr[i]   = (1.0)/outDeg[i];  
    }
    return;
}

void freeMem (graph* G)
{
    if (G->VI != NULL)
    {
        delete[] G->VI;
        G->VI = NULL;
    }
    if (G->EI != NULL)
    {
        delete[] G->EI;
        G->EI = NULL;
    }
    if (G->attr != NULL)
    {
        delete[] G->attr;
        G->attr = NULL;
    }
}
