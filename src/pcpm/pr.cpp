/**
 * Author: Kartik Lakhotia
 * Email id: klakhoti@usc.edu
 * Date: 27-Jul-2017
 *
 * This code implements work optimized propagation blocking with
 * transposed bin graph to reduce cache misses in scatter
 */


#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <omp.h>
#include <assert.h>
#include <immintrin.h>
#include "../../include/graph.hpp"
#ifndef SORT_HEADER
#include "../../include/sort.hpp"
#endif

using namespace std;


#define VERTEX_PARTITION
//#undef VERTEX_PARTITION

#define MAX_NEG 0x80000000
#define MAX_POS 0x7fffffff
#define MAX_UINT 0xffffffff

#define DEBUG
#undef DEBUG

#define DUMP
#undef DUMP

//////////////////////////////////////////
// level 2 debugging - asserts enabled
//////////////////////////////////////////
#define DEBUGL2
#undef DEBUGL2

//////////////////////////////////////////
// performance monitoring via PCM
//////////////////////////////////////////
#define PERF_MON
#undef PERF_MON

//////////////////////////////////////////
//parallel programming variables & types
//////////////////////////////////////////
typedef struct threadData
{
    int tid;
    int startVertex;
    int endVertex;
    graph* G;
}threadData;
int NUM_THREADS = 16;
int MAX_THREADS = 32;

//////////////////////////////////////////
//pagerank and binning variables
//////////////////////////////////////////
float dampingFactor = 0.15;
unsigned int MAX_ITER = 20;
unsigned int binWidth = (256*1024)/sizeof(float); //512kB
unsigned int binOffsetBits = (unsigned int)log2((float)binWidth); 
unsigned int NUM_BINS = 10000000/binWidth;

//////////////////////////////////////////
//define pre-processing functions
//////////////////////////////////////////
void partition(threadData*, graph*);
void computeOffsets(unsigned int**, unsigned int**, threadData*, graph*); 
void removeSameBinEdges(graph*); //deprecated
void transposePartition(graph*, graph*, threadData*, unsigned int*, unsigned int*);
template <class T> T** allocateBinMat (unsigned int);
template <class T> T*** allocateBinMatPtr (unsigned int);

//////////////////////////////////////////
//scatter-gather function definitions
//////////////////////////////////////////
void scatter(threadData*, float**, unsigned int**, unsigned int*, unsigned int*, unsigned int*); 
void gather(graph*, float*, unsigned int*, unsigned int*, unsigned int, unsigned int);
void scatterPull(threadData*, float*, float**);

//////////////////////////////////////////
////////////free the memory //////////////
//////////////////////////////////////////
template <class T> void freeMat (T**, unsigned int);


//////////////////////////////////////////
//main function
//////////////////////////////////////////
int main(int argc, char** argv)
{
    if (argc >= 3)
        NUM_THREADS = (unsigned int)atoi(argv[2]);
    if (argc >= 4)
        MAX_ITER = (unsigned int)atoi(argv[3]);
    if ((argc < 2) || (argc > 4))
    {
        printf("Usage : %s <filename> <numThreads(optional)> <#iterations(optional)> \n", argv[0]);
        exit(1);
    }
   
    omp_set_num_threads(NUM_THREADS);
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


#ifdef DEBUG
    printf("%d %d\n", G.numVertex, G.numEdges);
#endif    

    NUM_BINS = (G.numVertex-1)/binWidth + 1;

    //////////////////////////////////////////
    // output Degree array
    //////////////////////////////////////////
    unsigned int* outDeg = new unsigned int [G.numVertex]();
    findOutDeg(&G, outDeg);
#ifdef DEBUG
    printf("out degree computed\n");
#endif
    //////////////////////////////////////////
    // initialize page rank attribute to 1/degree
    //////////////////////////////////////////
    initGraphPullPR (&G, outDeg);
#ifdef DEBUG
    printf("page rank values initialized\n");
#endif 

    //////////////////////////////////////////
    //static work allocation to threads
    //equal no. of edges to all bins
    //////////////////////////////////////////
    threadData* TD = (threadData*) malloc (sizeof(threadData)*NUM_BINS);
    partition(TD, &G);

#ifdef DEBUGL2
    printf("static task allocation for scatter phase done\n");
    for (unsigned int i=0; i<NUM_BINS; i++)
    {
        printf("stats for partition %d:\n", i);
        printf("start_vertex=%d, end_vertex=%d, numEdges=%d\n", TD[i].startVertex, TD[i].endVertex, G.VI[TD[i].endVertex]-G.VI[TD[i].startVertex]);
    }
#endif 

    //////////////////////////////////////////
    //compute storage space required for each bin and
    //offsets for storage in bins for a thread
    //////////////////////////////////////////
    unsigned int** updateBinAddrOffset = allocateBinMat<unsigned int>(NUM_BINS+1);
    unsigned int** destIdBinAddrOffset = allocateBinMat<unsigned int>(NUM_BINS+1);

    struct timespec preStart, preEnd; 
    float preTime;
    if( clock_gettime(CLOCK_REALTIME, &preStart) == -1) { perror("clock gettime");}

    //////////////////////////////////////////
    //// transpose and compute offsets ///////
    //////////////////////////////////////////
    graph* GSort = new graph [NUM_BINS];

    #pragma omp parallel for schedule(dynamic) num_threads(MAX_THREADS)
    for (unsigned int i=0; i<NUM_BINS; i++)
        transposePartition(&G, &GSort[i], &TD[i], updateBinAddrOffset[i+1], destIdBinAddrOffset[i+1]);
    computeOffsets(updateBinAddrOffset, destIdBinAddrOffset, TD, &G);

    if( clock_gettime( CLOCK_REALTIME, &preEnd) == -1 ) { perror("clock gettime");}      
    preTime = (preEnd.tv_sec - preStart.tv_sec)+ (int)(preEnd.tv_nsec - preStart.tv_nsec)/1e9;
    printf("%s, preprocessing time - %lf\n", argv[1], preTime);

//////////////////////////////////////////
//////////////// BINNING ////////////////
//////////////////////////////////////////

    //////////////////////////////////////////
    //allocate space in bins to store  ///////
    //updates and destination vertex id's ////
    //////////////////////////////////////////
    float** updateBins = new float* [NUM_BINS];
    unsigned int** destIdBins = new unsigned int* [NUM_BINS];
    for (unsigned int i=0; i<NUM_BINS; i++)
    {
        updateBins[i] = new float [updateBinAddrOffset[NUM_BINS][i]];
        destIdBins[i] = new unsigned int [destIdBinAddrOffset[NUM_BINS][i]+1];
    }

    //////////////////////////////////////////
    ////individual pointers to starting //////
    //locations of bins for every thread /////
    //////////////////////////////////////////
    float*** indUpdateBins = allocateBinMatPtr<float>(NUM_BINS);
    unsigned int*** indDestIdBins = allocateBinMatPtr<unsigned int>(NUM_BINS);
    #pragma omp parallel for
    for (unsigned int i=0; i<NUM_BINS; i++)
    {
        for (unsigned int j=0; j<NUM_BINS; j++)
        {
            indUpdateBins[i][j] = updateBins[j] + updateBinAddrOffset[i][j];
            indDestIdBins[i][j] = destIdBins[j] + destIdBinAddrOffset[i][j];
        }
    }
    unsigned int** updateBinPointers = allocateBinMat<unsigned int>(NUM_BINS);
    unsigned int** destIdBinPointers = allocateBinMat<unsigned int>(NUM_BINS);

    //////////////////////////////////////////
    /////set last value of dest id -ve ///////
    ////// to indicate stop in gather phase///
    //////////////////////////////////////////
    #pragma omp parallel for
    for (unsigned int i=0; i<NUM_BINS; i++)
        destIdBins[i][destIdBinAddrOffset[NUM_BINS][i]] = MAX_NEG;

#ifdef DEBUG
    printf("thread offsets to individual bins recomputed\n");
#endif 

//////////////////////////////////////////
//////////// BINNING COMPLETE ////////////
//////////////////////////////////////////


    #pragma omp parallel for schedule(dynamic) num_threads(NUM_THREADS)
    for (unsigned int i=0; i<NUM_BINS; i++)
        scatter(&TD[i], indUpdateBins[i], indDestIdBins[i], updateBinPointers[i], destIdBinPointers[i], outDeg);

#ifdef DEBUG
    printf("first scatter finished\n");    
#endif

    #pragma omp parallel for schedule(dynamic) num_threads(NUM_THREADS)
    for (unsigned int i=0; i<NUM_BINS; i++)
        gather(&G, updateBins[i], destIdBins[i], outDeg, i, destIdBinAddrOffset[NUM_BINS][i]);


    //change the operand from original graph to PNG
    for (unsigned int i=0; i<NUM_BINS; i++)
        TD[i].G = &GSort[i];


    //free the memory.
    //original graph no more required
    delete[] G.VI;
    delete[] G.EI;
    G.VI = NULL;
    G.EI = NULL;

    
#ifdef DEBUG
    printf("graph with edges to bins constructed\n");
#endif

    //////////////////////////////////////////
    //compute MAX_ITER iterations of pagerank
    //////////////////////////////////////////
    struct timespec start, end, half; 
    float time;
	struct timespec scatterStart, scatterEnd;
	float scatterTime = 0.0;

    if( clock_gettime(CLOCK_REALTIME, &start) == -1) { perror("clock gettime");}
    unsigned int numIter = 0;


    while (numIter < MAX_ITER)
    {

    	if( clock_gettime(CLOCK_REALTIME, &scatterStart) == -1) { perror("clock gettime");}

        #pragma omp parallel for schedule(dynamic) num_threads(NUM_THREADS)
        for (unsigned int i=0; i<NUM_BINS; i++)
            scatterPull(&TD[i], G.attr, indUpdateBins[i]);

    	if( clock_gettime( CLOCK_REALTIME, &scatterEnd) == -1 ) { perror("clock gettime");}        
		scatterTime += (scatterEnd.tv_sec - scatterStart.tv_sec) + (float)(scatterEnd.tv_nsec - scatterStart.tv_nsec)/1e9; 
#ifdef DEBUGL2
        printf("scatter Time = %lf\n", scatterTime);
#endif
    
        #pragma omp parallel for schedule(dynamic) num_threads(NUM_THREADS)
        for (unsigned int i=0; i<NUM_BINS; i++)
            gather(&G, updateBins[i], destIdBins[i], outDeg, i, destIdBinAddrOffset[NUM_BINS][i]);

        numIter++;
    }

#ifdef DEBUG
    printf("pagerank computation finished\n");
#endif 

    if( clock_gettime( CLOCK_REALTIME, &end) == -1 ) { perror("clock gettime");}        
    time = (end.tv_sec - start.tv_sec)+ (int)(end.tv_nsec - start.tv_nsec)/1e9;
//    printf("%s, scatter Time = %lf, gather Time = %lf, time = %lf, ", argv[1], scatterTime, time-scatterTime, time);
    printf("%s, %lf\n", argv[1], time);



    
    //////////////////////////////////////////////
    //print top pagerank values for verification
    //////////////////////////////////////////////
#ifdef DUMP
    for (unsigned int i=0; i<G.numVertex; i++)
    {
        if (outDeg[i]>0)
            G.attr[i] = G.attr[i]*outDeg[i];
    }
    mergeSortWOkey<float>(G.attr, 0, G.numVertex-1);
    FILE* fdump = fopen("dumpPR.txt", "w");
    if (fdump == NULL)
    {
        fputs("file error\n", stderr);
        exit(1);
    }
    int printVertices = (G.numVertex>1000) ? 1000 : G.numVertex;
    for (int i=0; i<printVertices; i++)
        fprintf(fdump, "%lf\n", G.attr[i]);
    fclose(fdump);
#endif
    //////////////////////////////////////

    // free allocated memory//
    freeMem(&G);
    for (unsigned int i=0; i<NUM_BINS; i++)
        freeMem(&GSort[i]);
    delete[] GSort;
    delete[] outDeg;
    free(TD);
    freeMat<float>(updateBins, NUM_BINS);
    freeMat<unsigned int>(destIdBins, NUM_BINS);
    freeMat<unsigned int>(updateBinAddrOffset, NUM_BINS+1);
    freeMat<unsigned int>(destIdBinAddrOffset, NUM_BINS+1);
    freeMat<unsigned int>(updateBinPointers, NUM_BINS);
    freeMat<unsigned int>(destIdBinPointers, NUM_BINS);
    freeMat<float*>(indUpdateBins, NUM_BINS);
    freeMat<unsigned int*>(indDestIdBins, NUM_BINS);
    
    return 0;
}

//////////////////////////////////////////
//////// pre-processing functions ////////
//////////////////////////////////////////

//compute equal sized (vertex/edge) partitions
void partition(threadData* TD, graph* G)
{
    unsigned int numEdgesPerBin = (G->numEdges-1)/NUM_BINS + 1;
    unsigned int numVertexPerBin = (G->numVertex-1)/NUM_BINS + 1;
    int vcount = 0;
    for (int i=0; i<NUM_BINS; i++)
    {
        TD[i].tid = i;
        TD[i].G = G;
        TD[i].startVertex = vcount;
#ifndef VERTEX_PARTITION
        while(G->VI[vcount] < ((i+1)*numEdgesPerBin)) 
        {
            vcount++;
            if (vcount >= G->numVertex-1)
            {
                cout << vcount << " " << G->numVertex << endl;
                break;
            }
        }
#else
        vcount += numVertexPerBin;
#endif
        TD[i].endVertex = vcount;
        TD[i].endVertex = (TD[i].endVertex > G->numVertex) ? G->numVertex : TD[i].endVertex;
    }
}

//allocate BIN x BIN space for offsets
template <class T> T** allocateBinMat (unsigned int numRows)
{
    T** pointerMat;
    pointerMat = new T* [numRows];
    for (unsigned int i=0; i<numRows; i++)
        pointerMat[i] = new T [NUM_BINS]();
    return pointerMat;
}

//allocate BIN x BIN space for pointers
template <class T> T*** allocateBinMatPtr (unsigned int numRows)
{
    T*** pointerMat;
    pointerMat = new T** [numRows];
    for (unsigned int i=0; i<numRows; i++)
        pointerMat[i] = new T* [NUM_BINS];
    return pointerMat;
}


//compute individual partition offsets lock free reads/writes to memory
void computeOffsets(unsigned int** updateBinAddrOffset, unsigned int** destIdBinAddrOffset, threadData* TD, graph* G) 
{
    //accumulate the individual storage requirement
    //for different bins
    #pragma omp parallel for
    for (unsigned int i=0; i<NUM_BINS; i++)
    {
        for (unsigned int j=0; j<NUM_BINS; j++)
        {
            updateBinAddrOffset[j+1][i] += updateBinAddrOffset[j][i];
            destIdBinAddrOffset[j+1][i] += destIdBinAddrOffset[j][i];
        }
    }

#ifdef DEBUGL2
    printf("total edges incident on different bins:\n");
    for (unsigned int i=0; i<NUM_BINS; i++)
        printf("%d, %d\n", i, destIdBinAddrOffset[NUM_BINS][i]);
#endif

#ifdef DEBUG
    unsigned int redGraphSize = 0;
    for (unsigned int i=0; i<NUM_BINS; i++)
        redGraphSize += updateBinAddrOffset[NUM_BINS][i];
    printf("reduced graph size=%d\n", redGraphSize);
#endif 
}

//transpose the graph to sort on destination
void transposePartition(graph *G, graph* GSort, threadData* TD, unsigned int* updateBinAddrOffset, unsigned int* destIdBinAddrOffset)
{

    unsigned int currBin, prevBin;

    GSort->attr = NULL;

    GSort->numVertex = NUM_BINS;
    GSort->VI = new unsigned int [GSort->numVertex+1]();

    for (unsigned int i=TD->startVertex; i<TD->endVertex; i++)
    {
        prevBin = NUM_BINS+1;
        for (unsigned int j=G->VI[i]; j<G->VI[i+1]; j++)
        {
            currBin = (G->EI[j] >> binOffsetBits);
            destIdBinAddrOffset[currBin]++;
            if (currBin == prevBin)
                continue;
            GSort->VI[currBin+1]++;
            prevBin = currBin;
        }
    }

    for (unsigned int i=0; i<NUM_BINS; i++)
        updateBinAddrOffset[i] = GSort->VI[i+1];

    for (unsigned int i=0; i<GSort->numVertex; i++)
        GSort->VI[i+1] += GSort->VI[i];

    GSort->numEdges = GSort->VI[GSort->numVertex];
    GSort->EI = new unsigned int [GSort->numEdges];

    unsigned int* binOffset = new unsigned int [NUM_BINS]();
    for (unsigned int i=TD->startVertex; i<TD->endVertex; i++)
    {
        prevBin = NUM_BINS;
        for (unsigned int j=G->VI[i]; j<G->VI[i+1]; j++)
        {
            currBin = (G->EI[j] >> binOffsetBits);
            if (currBin == prevBin)
                continue; 
            
            GSort->EI[GSort->VI[currBin] + (binOffset[currBin]++)] = i;
            prevBin = currBin;
        }
    }

    delete[] binOffset;

}

//////////////////////////////////////////

//////////////////////////////////////////
//////// scatter gather functions ////////
//////////////////////////////////////////
void scatter(threadData* TD, float** updateBins, unsigned int** destIdBins, unsigned int* updateBinPointers, unsigned int* destIdBinPointers, unsigned int* outDeg)
{
    graph* G = TD->G;
    float prVal = 0.0;
    unsigned int destId = 0;
    unsigned int destBin = 0;
    unsigned int prevBin = 0;
    unsigned int loc = 0;
#ifdef DEBUGL2
    for (unsigned int i=0; i<NUM_BINS; i++)
    {
        assert(updateBinPointers[i] == 0);
        assert(destIdBinPointers[i] == 0);
    }
#endif
    for (unsigned int i=TD->startVertex; i<TD->endVertex; i++)
    {
        prevBin = NUM_BINS;
        prVal = G->attr[i];
        for (unsigned int j=G->VI[i]; j<G->VI[i+1]; j++)
        {
            destId = G->EI[j];
            destBin = (destId >> binOffsetBits);
            if (destBin!=prevBin)
            {
                updateBins[destBin][updateBinPointers[destBin]++] = prVal;
                destId |= MAX_NEG;
                prevBin = destBin;
            }
            destIdBins[destBin][destIdBinPointers[destBin]++] = destId;
        }
    }
    for (unsigned int i=0; i<NUM_BINS; i++)
    {
        destIdBinPointers[i] = 0;
        updateBinPointers[i] = 0;
    }
}


void scatterPull(threadData* TD, float* pr, float** updateBins)
{
    graph* G = TD->G;
    unsigned int pointer;
    unsigned int lastVal;
    for (unsigned int i=0; i<G->numVertex; i++)
    {
        pointer = 0;
        for (unsigned int j=G->VI[i]; j<G->VI[i+1]; j++)
        {
            updateBins[i][pointer++] = pr[G->EI[j]];
        }
    }  
}

void gather(graph* G, float* updateBin, unsigned int* destIdBin, unsigned int* outDeg, unsigned int binId, unsigned int binSize)
{
    unsigned int startVertex = (binId << binOffsetBits);
    unsigned int endVertex = ((binId+1) << binOffsetBits);
    endVertex = (endVertex > G->numVertex) ? G->numVertex : endVertex;
    for (unsigned int i=startVertex; i<endVertex; i++)
        G->attr[i] = 0;
    unsigned int destBinPointer = 0;
    unsigned int updateBinPointer = MAX_UINT;
    float prValRecord = 0.0;
    unsigned int destId = 0;
    for (unsigned int i=0; i<binSize; i++)
    {
#ifdef DEBUGL2
        assert((destIdBin[i] & MAX_POS) < endVertex);
        assert((destIdBin[i] & MAX_POS) >= startVertex);
#endif
        destId = destIdBin[i];
        updateBinPointer += (destId >> 31);
        G->attr[destId & MAX_POS] += updateBin[updateBinPointer]; 
    }

    for (unsigned int i=startVertex; i<endVertex; i++)
    {
        if (outDeg[i] > 0)
            G->attr[i] = (dampingFactor + (1-dampingFactor)*G->attr[i])/(outDeg[i]);
    }


}

//////////////////////////////////////////
////////////free the memory //////////////
//////////////////////////////////////////
template <class T> void freeMat (T** mat, unsigned int numRows)
{
    for (unsigned int i=0; i<numRows; i++)
        delete[] mat[i];
    delete[] mat; 
}
