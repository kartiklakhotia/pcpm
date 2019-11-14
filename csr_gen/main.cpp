#include <iostream>
#include <vector>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <unordered_map>
#include "graph.h"

#define DEBUG 
#undef DEBUG

bool weighted = false;
bool createReverse = false;
bool undirected = false;

using namespace std;

int main(int argc, char** argv)
{

    unsigned int wtValue = 0;
    unsigned int reverseFileIndex = 0;

    if (argc < 3)
    {
        printf("Usage : %s <inputFile1> <outputFile> -w <weight_type> -r <graphTransposeOutputFilename>\n", argv[0]);
        exit(1);
    }
    for (int i=1; i<argc; i++)
    {
        if (strcmp(argv[i], "-w")==0)
        {
            weighted = true;
            if (i+1 < argc)
                wtValue = (unsigned int)atoi(argv[i+1]);
            else
            {
                printf("specify whether to read weights from file or assign all 1's\n");
                exit(1);
            }
        }
        if (strcmp(argv[i], "-r")==0)
        {
            createReverse = true;
            reverseFileIndex = i+1;
        }
        if (strcmp(argv[i], "-u")==0)
            undirected = true;
    }

    // edge list
    std::vector<intV> src;
    std::vector<intV> dst;
    std::vector<unsigned int> ew;

    intE numEdgesRead = 0;

    FILE* fp = fopen (argv[1], "r");
    if (fp == NULL)
    {
        fputs("file error", stderr);
        return -1;
    }

    intV srcVal, dstVal;
    intV numVertex = 0;
    unsigned int edgeWeight = 0;
    while(!feof(fp))
    {
        if(fscanf(fp, "%d", &srcVal) <= 0)
            break;
        fscanf(fp, "%d", &dstVal);
        if (weighted)
        {
            if (wtValue == 0)
                fscanf(fp, "%d", &edgeWeight);
            else if (wtValue == 1)
                edgeWeight = 1;
            else
                edgeWeight = rand() % 10 + 1;
        }
       
        numVertex = (srcVal > numVertex) ? srcVal : numVertex;
        numVertex = (dstVal > numVertex) ? dstVal : numVertex;
        if (srcVal != dstVal)
        {
            src.push_back(srcVal);
            dst.push_back(dstVal);
            ew.push_back(edgeWeight);
            numEdgesRead++;
            if (undirected)
            {
                src.push_back(dstVal);
                dst.push_back(srcVal);
                ew.push_back(edgeWeight);
                numEdgesRead++;
            }
        }
    }
    fclose(fp);

    numVertex++;

    intV* inDeg = new intV [numVertex]();

    
    for (intE i=0; i<numEdgesRead; i++)
        inDeg[src[i]]++;


    graph G1;
    G1.weighted = weighted;
    G1.numVertex = numVertex;
    G1.numEdges = numEdgesRead;
    G1.VI = new intE [numVertex+1]();
    G1.EI = new intV [numEdgesRead]();
    if (weighted)
        G1.EW = new unsigned int [numEdgesRead]();

    for (intV i=1; i<=numVertex; i++)
        G1.VI[i] = G1.VI[i-1] + inDeg[i-1];

    for (intV i=0; i<numVertex; i++)
        inDeg[i] = 0;

    for (intE i=0; i<numEdgesRead; i++)
    {
        G1.EI[G1.VI[src[i]] + inDeg[src[i]]] = dst[i];
        if (weighted)
            G1.EW[G1.VI[src[i]] + inDeg[src[i]]] = ew[i];
        inDeg[src[i]]++;
    }

    sortEdges(&G1);
    write_csr(argv[2], &G1); 

    if (createReverse)
    {
        graph G2;
        createReverseCSR(&G1, &G2);
        write_csr(argv[reverseFileIndex], &G2); 
        freeMem(&G2);
    }


    delete[] inDeg;
    freeMem(&G1);

#ifdef DEBUG
    read_csr (argv[2], &G1);
    printGraph(&G1);
    freeMem(&G1);
    if (createReverse)
    {
        read_csr (argv[reverseFileIndex], &G1);
        printGraph(&G1);
        freeMem(&G1);
    }
#endif

    return 0;
}

