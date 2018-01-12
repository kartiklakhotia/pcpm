/**
 * Author: Kartik Lakhotia
 * Email id: klakhoti@usc.edu
 * Date: 22-Sep-2017
 *
 * This code implements work optimized propagation blocking with
 * transposed bin graph to reduce cache misses in scatter
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <vector>
#include <queue>
#include <string>
#include <algorithm>
#include <omp.h>
#include <assert.h>
#include <immintrin.h>
#ifndef GRAPH_HEADER
#include "../../include/graph.hpp"
#endif
#ifndef SORT_HEADER
#include "../../include/sort.hpp"
#endif
#include "../../include/partFunc.hpp"

using namespace std;


unsigned int order = 3; // 0 -> original, 1 -> ascending degree, 2 -> descending degree, 3 -> random
unsigned int PART_SIZE   = (256*1024)/sizeof(float); //256kB
unsigned int TARGET_COMP = 5;

bool DUMP_MAP = false;


int main(int argc, char** argv)
{
    if ((argc<3) || (argc>5))
    {
        printf("Usage : %s <inputFilename> <outputFileName> <traversalOrder(Optional)> <dumpNewLabels(optional)>\n", argv[0]);
        exit(1);
    }
    if (argc > 3)
        order = (unsigned int)stoi(argv[3]);
    if (argc > 4)
        DUMP_MAP = (stoi(argv[4])==0) ? false : true; 


    // graph object
    graph G;

    //////////////////////////////////////////
    // read csr file
    //////////////////////////////////////////
    if (read_csr(argv[1], &G)==-1)
    {
        printf("couldn't read %s\n", argv[1]);
        exit(1);
    }
    G.attr = NULL;


    unsigned int NUM_PARTS = (G.numVertex-1)/PART_SIZE + 1;
    unsigned int LAST_PART = NUM_PARTS - 1;
    unsigned int padding = NUM_PARTS*PART_SIZE - G.numVertex;
 
    //////////////////////////////////////////
    // output Degree array
    //////////////////////////////////////////
    unsigned int* outDeg = new unsigned int [G.numVertex]();
    findOutDeg(&G, outDeg);

    unsigned int* nodeId = new unsigned int [G.numVertex]();
    for (unsigned int i=0; i<G.numVertex; i++)
        nodeId[i] = i;


    // obtain order of traversal of nodes based on specification //
    if (order == 0) {;}
    else if (order == 1) 
        mergeSortAscending<unsigned int, unsigned int>(outDeg, nodeId, 0, G.numVertex-1);
    else if (order == 2)
        mergeSort<unsigned int, unsigned int>(outDeg, nodeId, 0, G.numVertex-1);
    else if (order == 3)
        randomize(nodeId, outDeg, G.numVertex);



    unsigned int degThresh = (unsigned int)(logProd<unsigned int>(TARGET_COMP)*((double)NUM_PARTS));

    /////////////////////////////////////////////
    //Each node is associated with partition ID//
    //and a sequence ID indicative of when it  //
    //was assigned that partition       /////////
    /////////////////////////////////////////////

    //array to store partition assignments
    unsigned int* partId = new unsigned int [G.numVertex];
    for (unsigned int i=0; i<G.numVertex; i++)
        partId[i] = NUM_PARTS; //to indicate that no partition has been assigned yet

    //order in which nodes are assigned to the partition.
    unsigned int* seqId = new unsigned int [G.numVertex]();
     

    //node to partition edge weight for assigned neighbors
    unsigned int* n2pWt = new unsigned int [NUM_PARTS]();

    //list of partitions that are neighbors(and degree)
    unsigned int* np = new unsigned int [NUM_PARTS]();
    unsigned int dp = 0;


    /////////////////////////////////////////////////
    /// partition heap indicating number of nodes ///
    /////////////////////////////////////////////////
    heap* partCont = createHeap(NUM_PARTS);
    initializeHeap(partCont, G.numVertex);
    partCont->elems[partCont->size-1]->numNodes = padding;

    unsigned int minPart, minPartNodes;



    /////////////////////////////////////////////////
    ///////// compute partition assignment //////////    
    /////////////////////////////////////////////////
    for (unsigned int i=0; i<G.numVertex; i++)
    {
        // create list of neighboring partitions
        if ((outDeg[i] > 1) && (outDeg[i] < degThresh)) //node with degree 1 will not compress and with >degThresh will exhibit target compression always
        {
            dp = createNeighborPartList (&G, nodeId[i], partId, NUM_PARTS, n2pWt, np);

            
            //find least filled partition the unassigned neighbors should belong to
            minPart = leastFilledNeighborPart (NUM_PARTS, PART_SIZE, np, dp, partCont, n2pWt);
            if (minPart == NUM_PARTS) //no partition was selected because no neighbors assigned or all neighbor partitions full
                minPart = getMinPartId(partCont);

            minPartNodes = getNumNodes(partCont, minPart);

            //assign the remaining nodes to previously computed partition
            assignNodes(&G, nodeId[i], partId, NUM_PARTS, PART_SIZE, seqId, n2pWt, minPart, minPartNodes, partCont, np, dp);
        }
    }



    //just in case some nodes were left unassigned
    //assign them to least filled partitions
    minPart = 0;
    minPartNodes = getNumNodes(partCont, minPart);
    for (unsigned int i=0; i<G.numVertex; i++)
    {
        if (partId[i]==NUM_PARTS)
        {
            while (minPartNodes >= PART_SIZE)
            {
                minPart++;
                minPartNodes = getNumNodes(partCont, minPart);
            }
            partId[i] = minPart;
            seqId[i] = minPartNodes++;
            updateKey(partCont, minPart, minPartNodes);
        }
    }



    ////recompute number of edges////
    unsigned int numCompressedEdges = 0;
    for (unsigned int i=0; i<G.numVertex; i++)
    {
        numCompressedEdges += createNeighborPartList(&G, i, partId, NUM_PARTS, n2pWt, np); 
        for (unsigned int j=G.VI[i]; j<G.VI[i+1]; j++)
            n2pWt[partId[G.EI[j]]] = 0;
    }

    printf("number of compressed edges is %d\n", numCompressedEdges);



    /////////////////////////////////////
    ///////////node labelling////////////
    /////////////////////////////////////

    //compute prefix sum to find the label
    //of first node in each partition
    unsigned int* prefix = new unsigned int[NUM_PARTS]();
    for (unsigned int i=1; i<NUM_PARTS; i++)
    { 
        prefix[i] = prefix[i-1] + getNumNodes(partCont, i-1);
    }

    // heap won't be used anymore //
    freeHeapMem(partCont);

    //compute new node labels
    unsigned int* newNodeId = new unsigned int[G.numVertex]; // oldId->newId mapping
    for (unsigned int i=0; i<G.numVertex; i++)
    {
        newNodeId[i] = prefix[partId[i]] + seqId[i];
        if (partId[i] == LAST_PART)
            newNodeId[i] = newNodeId[i] - padding;
    }

    if (DUMP_MAP)
        dumpNewOrder (newNodeId, G.numVertex);        

    //construct new graph//
    graph Gnew;
    applyNewLabeling(&G, &Gnew, newNodeId);


    //free allocated memory//
    //write the graph in new file//
    write_csr(argv[argc-1], &Gnew);

    freeMem(&G);
    delete[] outDeg;
    delete[] nodeId;
    delete[] n2pWt;
    delete[] np;
    delete[] newNodeId;
    delete[] partId;
    delete[] seqId;
    freeMem(&Gnew);

    return 0;
}
