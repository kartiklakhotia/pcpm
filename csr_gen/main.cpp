#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <unordered_map>
#include "graph.h"

#define DEBUG 
#undef DEBUG

using namespace std;



int main(int argc, char** argv)
{

    if (argc != 3)
    {
        printf("Usage : %s <inputFile1> <outputFile>\n", argv[0]);
        exit(1);
    }

    // edge list
    std::vector<unsigned int> src;
    std::vector<unsigned int> dst;

    unsigned int numEdgesRead = 0;

    FILE* fp = fopen (argv[1], "r");
    if (fp == NULL)
    {
        fputs("file error", stderr);
        return -1;
    }

    unsigned int srcVal, dstVal;
    unsigned int numVertex = 0;
//    unsigned int edgeWeight = 0;
    while(!feof(fp))
    {
        if(fscanf(fp, "%d", &srcVal) <= 0)
            break;
        fscanf(fp, "%d", &dstVal);
//        fscanf(fp, "%d", &edgeWeight);
        numVertex = (srcVal > numVertex) ? srcVal : numVertex;
        numVertex = (dstVal > numVertex) ? dstVal : numVertex;
        if (srcVal != dstVal)
        {
            src.push_back(srcVal);
            dst.push_back(dstVal);
            numEdgesRead++;
        }
    }
    fclose(fp);

    numVertex++;

    unsigned int* inDeg = new unsigned int [numVertex]();

    
    for (unsigned int i=0; i<numEdgesRead; i++)
        inDeg[src[i]]++;


    graph G1;
    G1.numVertex = numVertex;
    G1.numEdges = numEdgesRead;
    G1.VI = new unsigned int [numVertex+1]();
    G1.EI = new unsigned int [numEdgesRead]();

    for (unsigned int i=1; i<=numVertex; i++)
        G1.VI[i] = G1.VI[i-1] + inDeg[i-1];

    for (unsigned int i=0; i<numVertex; i++)
        inDeg[i] = 0;

    for (unsigned int i=0; i<numEdgesRead; i++)
    {
        G1.EI[G1.VI[src[i]] + inDeg[src[i]]] = dst[i];
        inDeg[src[i]]++;
    }

    sortEdges(&G1);
    write_csr(argv[2], &G1); 


    delete[] inDeg;
    freeMem(&G1);

#ifdef DEBUG
    read_csr (argv[2], &G1);
    printGraph(&G1);
    freeMem(&G1);
#endif

    return 0;
}

