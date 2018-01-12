#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#ifndef GRAPH_HEADER
#include "graph.hpp"
#endif


// datatype for partition min-heap
// min-heap used to find least filled partition
typedef struct heapNode
{
    unsigned int partId;
    unsigned int numNodes;
} heapNode;

//heap structure and related functions to find
//least filled partition
typedef struct heap
{
    unsigned int size;
    heapNode** elems;         
    unsigned int* pos; //indicate position of a partId in elems array
    unsigned int maxNodes; //to indicate infinity. No partition can have more than maxNodes
} heap;

heap* createHeap(unsigned int size)
{
    heap* heapInstance = new heap;
    heapInstance->size = size + 1;
    heapInstance->pos = new unsigned int [size]();
    heapInstance->elems = new heapNode* [heapInstance->size];
    for (unsigned int i=0; i<heapInstance->size; i++)
        heapInstance->elems[i] = new heapNode;
    return heapInstance;
}

void initializeHeap(heap* heapInstance, unsigned int maxNodes)
{
    heapInstance->maxNodes = maxNodes;

    for (unsigned int i=1; i<heapInstance->size; i++)
    {
        heapInstance->pos[i-1] = i;
        heapInstance->elems[i]->partId = i-1;
        heapInstance->elems[i]->numNodes = 0;
    }
}

//swap 2 nodes at posA, posB in elems array
void swapNodes(heap* heapInstance, unsigned int posA, unsigned int posB)
{
    unsigned int partAId = heapInstance->elems[posA]->partId;
    unsigned int partBId = heapInstance->elems[posB]->partId; 
    heapNode* temp = heapInstance->elems[posA];
    heapInstance->elems[posA] = heapInstance->elems[posB];
    heapInstance->elems[posB] = temp;
    heapInstance->pos[partAId] = posB;
    heapInstance->pos[partBId] = posA;
}

unsigned int getNumNodes (heap* heapInstance, unsigned int partId)
{
    return heapInstance->elems[heapInstance->pos[partId]]->numNodes;
}

unsigned int getMinPartId (heap* heapInstance)
{
    return heapInstance->elems[1]->partId;
}

void updateKey(heap* heapInstance,  unsigned int partId, unsigned int newNodes)
{
    unsigned int partPos = heapInstance->pos[partId];
    heapInstance->elems[partPos]->numNodes = newNodes;
    unsigned int childPosL, childPosR, childNodesL, childNodesR;
    unsigned int childPosMin, childNodesMin;

    //swap until reach bottom or when both children are greater
    while(1)
    {
        partPos = heapInstance->pos[partId];
        childPosL = partPos*2;
        childPosR = partPos*2+1;
        if ((childPosL >= heapInstance->size) && (childPosR >= heapInstance->size))
            return;
        childNodesL = (childPosL < heapInstance->size) ? heapInstance->elems[childPosL]->numNodes : heapInstance->maxNodes;
        childNodesR = (childPosR < heapInstance->size) ? heapInstance->elems[childPosR]->numNodes : heapInstance->maxNodes;
        if (childNodesL <= childNodesR)
        {
            childPosMin = childPosL;
            childNodesMin = childNodesL;
        }
        else
        {
            childPosMin = childPosR;
            childNodesMin = childNodesR;
        }
        if (childNodesMin >= newNodes)
            return;
        swapNodes(heapInstance, partPos, childPosMin);
    }
} 

void printHeap(heap* heapInstance)
{
    printf("printing heap in (partId, numNodes) format\n");
    for (unsigned int i=1; i<heapInstance->size; i=i*2)
    {
        for (unsigned j=i; j<i*2; j++)
            printf("(%d, %d) ", heapInstance->elems[j]->partId, heapInstance->elems[j]->numNodes);
        printf("\n"); 
    }
}

heapNode* pop(heap* heapInstance)
{
    heapNode* retNode = heapInstance->elems[1];
    unsigned int lastPos = heapInstance->size-1;
    unsigned int lastNodes = heapInstance->elems[lastPos]->numNodes;
    unsigned int lastPartId = heapInstance->elems[lastPos]->partId;
    swapNodes(heapInstance, 1, lastPos);
    heapInstance->size = heapInstance->size-1;
    updateKey(heapInstance, lastPartId, lastNodes);  
}

void freeHeapMem(heap* heapInstance)
{
    for (unsigned int i=0; i<heapInstance->size; i++)
        delete[] heapInstance->elems[i];
    delete[] heapInstance->pos;
    delete[] heapInstance->elems;
    delete[] heapInstance;
}


// returns r*W(-r*exp(-r)) where r -> target compression ratio
// and W -> lambert function.
// Refer to partitioning.docx for more details
template <typename T>
double logProd (T r)
{
    double thresh = 0;
    if (r<=1)
        thresh = 0;
    else if (r<=2)
        thresh = 1.59;
    else if (r<=3)
        thresh = 2.82;
    else if (r<=4)
        thresh = 3.92;
    else if (r<=5)
        thresh = 4.96;
    else if (r<=6)
        thresh = 5.98;
    else
        thresh = r;
    return thresh; 
}

void randomize (unsigned int* nodeId, unsigned int* outDeg, unsigned int numVertex)
{
    std::vector<unsigned int> v;
    v.reserve(numVertex);
    for (unsigned int i=0; i<numVertex; i++)
        v.push_back(nodeId[i]);
    std::random_shuffle(v.begin(), v.end());
    for (unsigned int i=0; i<numVertex; i++)
    {
        nodeId[i] = v[i];
        v[i] = outDeg[i];
    }
    for (unsigned int i=0; i<numVertex; i++)
        outDeg[i] = v[nodeId[i]];
}

//create list of neighboring partitions
unsigned int createNeighborPartList (graph* G, unsigned int nodeId, unsigned int* partId, unsigned int numParts, unsigned int* n2pWt, unsigned int* list)
{
    unsigned int neighbor, neighborPart;
    unsigned int listSize = 0;
    for (unsigned int i=G->VI[nodeId]; i<G->VI[nodeId+1]; i++)
    {
        unsigned int neighbor = G->EI[i];
        unsigned int neighborPart = partId[neighbor];
        if (neighborPart < numParts) //neighbor has been assigned to a partition already
        {
            if (n2pWt[neighborPart]==0) //partition hasn't been added to list yet
                list[listSize++]=neighborPart;
            n2pWt[neighborPart]++;
        }
    }
    return listSize;
}


//find which partition the remaining nodes should belong to
unsigned int leastFilledNeighborPart (unsigned int numParts, unsigned int partSize, unsigned int* list, unsigned int listSize, heap* heapInstance, unsigned int* weight)
{
    unsigned int minPart = numParts;
    unsigned int minPartNodes = heapInstance->maxNodes;
    unsigned int maxWt = 0;
    unsigned int neighborPart, neighborPartNodes;
    for (unsigned int i=0; i<listSize; i++)
    {
        neighborPart = list[i];
        neighborPartNodes = getNumNodes(heapInstance, neighborPart);
//        if ((neighborPartNodes < minPartNodes) && (neighborPartNodes < partSize))
        if ((weight[neighborPart] > maxWt) && (neighborPartNodes < partSize))
        {
            maxWt = weight[neighborPart];
            minPartNodes = neighborPartNodes;
            minPart = neighborPart;
        }
    }
    return minPart;
}

//assign the remaining nodes to partitions. return the number of partitions that were used
void assignNodes(graph*G, unsigned int nodeId, unsigned int* partId, unsigned int numParts, unsigned int partSize, unsigned int* seqId, unsigned int* n2pWt, unsigned int minPart, unsigned int minPartNodes, heap* heapInstance, unsigned int* list, unsigned int listSize)
{
    unsigned int neighbor, neighborPart;
    for (unsigned int i=G->VI[nodeId]; i<G->VI[nodeId+1]; i++)
    {
        neighbor = G->EI[i];
        neighborPart = partId[neighbor];
        if (neighborPart >= numParts)
        {
            partId[neighbor] = minPart;
            seqId[neighbor] = minPartNodes++;
            if (minPartNodes >= partSize) //if the partition gets full
            {
                updateKey(heapInstance, minPart, minPartNodes);
                minPart = leastFilledNeighborPart (numParts, partSize, list, listSize, heapInstance, n2pWt);
                if (minPart == numParts)
                    minPart = getMinPartId(heapInstance);
                minPartNodes = getNumNodes(heapInstance, minPart);
            }
        }
    }
    
    //update number of nodes of minPart in heap
    updateKey(heapInstance, minPart, minPartNodes);


    //reset the weights to 0
    for (unsigned int i=G->VI[nodeId]; i<G->VI[nodeId+1]; i++)
        n2pWt[partId[G->EI[i]]] = 0; 

}

void applyNewLabeling(graph* Gold, graph* Gnew, unsigned int* newNodeId)
{
    unsigned int* oldNodeId = new unsigned int[Gold->numVertex]; //newId->oldId mapping
    for (unsigned int i=0; i<Gold->numVertex; i++)
        oldNodeId[newNodeId[i]] = i;

    Gnew->numVertex = Gold->numVertex;
    Gnew->numEdges = Gold->numEdges; 
    Gnew->attr = NULL;
    Gnew->VI = new unsigned int [Gnew->numVertex+1]();
    Gnew->EI = new unsigned int [Gnew->numEdges](); 
    unsigned int base, offset;
    for (unsigned int i=0; i<Gnew->numVertex; i++)
    {
        Gnew->VI[i+1] = Gnew->VI[i] + (Gold->VI[oldNodeId[i]+1] - Gold->VI[oldNodeId[i]]);
        base = Gold->VI[oldNodeId[i]]; 
        offset = 0;
        for (unsigned int j=Gnew->VI[i]; j<Gnew->VI[i+1]; j++)
            Gnew->EI[j] = newNodeId[Gold->EI[base + offset++]];
    }
    sortEdges(Gnew);



    delete[] oldNodeId;
} 

void dumpNewOrder (unsigned int* newNodeId, unsigned int size)
{
    FILE* fp = fopen("order.txt", "w");
    for (unsigned int i=0; i<size; i++)
        fprintf(fp, "%d\n", newNodeId[i]);
}
