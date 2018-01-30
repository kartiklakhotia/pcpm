/**
 * Author: Kartik Lakhotia
 * Email id: klakhoti@usc.edu
 * Date: 5-Jun-2017
 *
 * This code implements propagation blocking to compute parallel pagerank
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
#include "../../include/utils.hpp"
#ifndef SORT_HEADER
#include "../../include/sort.hpp"
#endif

using namespace std;


#define DEBUG
#undef DEBUG

#define DUMP
//#undef DUMP

#define DEBUGL2
#undef DEBUGL2

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

//////////////////////////////////////////
//pagerank and binning variables
//////////////////////////////////////////
float dampingFactor = 0.15;
unsigned int binWidth = (512*1024)/sizeof(float); //512kB
unsigned int binOffsetBits = (unsigned int)log2((float)binWidth); 
unsigned int bufferSize = (2*64)/(sizeof(float)); 
unsigned int streamSize = 256/(8*sizeof(float));
unsigned int numStreams = bufferSize/streamSize;
unsigned int MAX_ITER = 10;
unsigned int NUM_BINS = 10000000/binWidth;

//////////////////////////////////////////
//scatter-gather function definitions
//////////////////////////////////////////
void scatter(threadData*, float**, unsigned int**, unsigned int*, unsigned int*); 
void scatterDPB(threadData*, float**,  float**, unsigned int*, unsigned int*, unsigned int*);
//void gather(graph*, float*, unsigned int*, unsigned int, unsigned int);

//////////////////////////////////////////
//// for separate bins to each thread ////
//////////////////////////////////////////
void gather(graph*, float***, unsigned int***, unsigned int, unsigned int**); 


//////////////////////////////////////////
////////// main function /////////////////
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

//    sortEdges(&G);
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
    initGraph (&G, outDeg);
#ifdef DEBUG
    printf("page rank values initialized\n");
#endif 

    //////////////////////////////////////////
    //static work allocation to threads
    //equal no. of edges to all bins
    //////////////////////////////////////////
    threadData* TD = (threadData*) malloc (sizeof(threadData)*NUM_THREADS);
    unsigned int numEdgesPerThread = (G.numEdges-1)/NUM_THREADS + 1;
    int vcount = -1;
    for (int i=0; i<NUM_THREADS; i++)
    {
        TD[i].tid = i;
        TD[i].G = &G;
        TD[i].startVertex = vcount + 1;
        vcount = vcount + 1;
        while(G.VI[vcount] < ((i+1)*numEdgesPerThread)) 
        {
            vcount++;
            if (vcount == G.numVertex-1)
                break;
        }
        TD[i].endVertex = vcount+1;
    }

#ifdef DEBUG
    printf("static task allocation for scatter phase done\n");
    for (unsigned int i=0; i<NUM_THREADS; i++)
    {
        printf("stats for thread %d:\n", i);
        printf("start_vertex=%d, end_vertex=%d, numEdges=%d\n", TD[i].startVertex, TD[i].endVertex, G.VI[TD[i].endVertex]-G.VI[TD[i].startVertex]);
    }
#endif 

    struct timespec preStart, preEnd; 
    float preTime;

    if( clock_gettime(CLOCK_REALTIME, &preStart) == -1) { perror("clock gettime");}
    //////////////////////////////////////////
    //compute storage space required for each bin and
    //offsets for storage in bins for a thread
    //////////////////////////////////////////
    NUM_BINS = (G.numVertex-1)/binWidth + 1;
    unsigned int** binAddrOffset = new unsigned int* [NUM_THREADS+1];
    for (unsigned int i=0; i<=NUM_THREADS; i++)
        binAddrOffset[i] = new unsigned int [NUM_BINS]();

    //first compute individual space of each thread in parallel    
    #pragma omp parallel for
    for (unsigned int i=0; i<NUM_THREADS; i++)
    {
        for (unsigned int j=G.VI[TD[i].startVertex]; j<G.VI[TD[i].endVertex]; j++)
        {
            //////////////////////////////////////////
            /// can fix the destination vertex ids ///
            //////////////////////////////////////////
            binAddrOffset[i+1][(G.EI[j] >> binOffsetBits)]++;    
        }
    }

    if( clock_gettime( CLOCK_REALTIME, &preEnd) == -1 ) { perror("clock gettime");}      
    preTime = (preEnd.tv_sec - preStart.tv_sec)+ (int)(preEnd.tv_nsec - preStart.tv_nsec)/1e9;
    printf("pre-processing time = %lf\n", preTime);
#ifdef DEBUG
    printf("bin space per thread computed\n");
    printf("number of bins is %d\n", NUM_BINS);
    for (unsigned int i=0; i<NUM_THREADS; i++)
    {
        unsigned int space = 0;
        for (unsigned int j=0; j<NUM_BINS; j++)
            space += binAddrOffset[i+1][j];
        printf("space alloted to thread %d is %d\n", i, space);
    }
#endif 
#ifdef DEBUG
    printf("bins' storage space allocated\n");
#endif 

    //////////////////////////////////////////
    //individual pointers to starting 
    //locations of bins for every thread
    //////////////////////////////////////////
    float*** indUpdateBins = new float** [NUM_THREADS];
    float*** buffers = new float** [NUM_THREADS];
    unsigned int*** indDestIdBins = new unsigned int** [NUM_THREADS];
    unsigned int** binPointers = new unsigned int* [NUM_THREADS];
    unsigned int** bufferPointers = new unsigned int* [NUM_THREADS];
    for (unsigned int i=0; i<NUM_THREADS; i++)
    {
        indUpdateBins[i] = new float* [NUM_BINS];
        indDestIdBins[i] = new unsigned int* [NUM_BINS];
        buffers[i] = new float* [NUM_BINS];
        for (unsigned int j=0; j<NUM_BINS; j++)
        {
            //////////////////////////////////////////////////////
            /// separate bins if bins are to be memory aligned ///
            //////////////////////////////////////////////////////
            indUpdateBins[i][j] = (float*) _mm_malloc(binAddrOffset[i+1][j]*sizeof(float), 32); 
            indDestIdBins[i][j] = new unsigned int [binAddrOffset[i+1][j]];
            buffers[i][j] = (float*) _mm_malloc(bufferSize*sizeof(float), 32);
        }
        binPointers[i] = new unsigned int[NUM_BINS]();
        bufferPointers[i] = new unsigned int [NUM_BINS]();
    }
#ifdef DEBUG
    printf("thread offsets to individual bins computed and stored\n");
#endif 
    


    #pragma omp parallel for schedule(static,1) num_threads(NUM_THREADS)
    for (unsigned int i=0; i<NUM_THREADS; i++)
        scatter(&TD[i], indUpdateBins[i], indDestIdBins[i], binPointers[i], outDeg);
    #pragma omp parallel for schedule(dynamic) num_threads(NUM_THREADS)
    for (unsigned int i=0; i<NUM_BINS; i++)
      gather(&G, indUpdateBins, indDestIdBins, i, binAddrOffset);



    //////////////////////////////////////////
    //compute MAX_ITER iterations of pagerank
    //////////////////////////////////////////
    struct timespec start, end; 
    float time;
	struct timespec scatterStart, scatterEnd;
	float scatterTime = 0.0;

    if( clock_gettime(CLOCK_REALTIME, &start) == -1) { perror("clock gettime");}
    unsigned int numIter = 0;


    while (numIter < MAX_ITER)
    {
    	if( clock_gettime(CLOCK_REALTIME, &scatterStart) == -1) { perror("clock gettime");}
        #pragma omp parallel for schedule(static,1) num_threads(NUM_THREADS)
        for (unsigned int i=0; i<NUM_THREADS; i++)
            scatterDPB(&TD[i], indUpdateBins[i], buffers[i], binPointers[i], bufferPointers[i], outDeg);


		
    	if( clock_gettime( CLOCK_REALTIME, &scatterEnd) == -1 ) { perror("clock gettime");}        
		scatterTime += (scatterEnd.tv_sec - scatterStart.tv_sec) + (float)(scatterEnd.tv_nsec - scatterStart.tv_nsec)/1e9; 

        #pragma omp parallel for schedule(dynamic) num_threads(NUM_THREADS)
        for (unsigned int i=0; i<NUM_BINS; i++)
          gather(&G, indUpdateBins, indDestIdBins, i, binAddrOffset);

        numIter++;
    }


#ifdef DEBUG
    printf("pagerank computation finished\n");
#endif 

    if( clock_gettime( CLOCK_REALTIME, &end) == -1 ) { perror("clock gettime");}        
    time = (end.tv_sec - start.tv_sec)+ (int)(end.tv_nsec - start.tv_nsec)/1e9;

    printf("%s, %lf\n", argv[1], time);

    //////////////////////////////////////////////
    //print top pagerank values for verification
    //////////////////////////////////////////////
#ifdef DUMP
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
    free(TD);
    for (unsigned int i=0; i<NUM_BINS; i++)
    {
        for (unsigned int j=0; j<NUM_THREADS; j++) 
        {
            free(indUpdateBins[j][i]);
            delete[] indDestIdBins[j][i];
            free(buffers[j][i]);
        }
    } 
    for (unsigned int i=0; i<NUM_THREADS; i++) 
    {
        delete[] binAddrOffset[i];
        delete[] indUpdateBins[i];
        delete[] indDestIdBins[i];
        delete[] binPointers[i];
        delete[] buffers[i];
    }

    delete[] binPointers;
    delete[] indUpdateBins;
    delete[] indDestIdBins;
    delete[] binAddrOffset;
    delete[] outDeg;
    delete[] buffers;
    
    

    return 0;
}

void scatter(threadData* TD, float** updateBins, unsigned int** destIdBins, unsigned int* binPointers, unsigned int* outDeg)
{
    graph* G = TD->G;
    float prVal = 0.0;
    unsigned int destId = 0;
    unsigned int destBin = 0;
    unsigned int loc = 0;
#ifdef DEBUGL2
    for (unsigned int i=0; i<NUM_BINS; i++)
        assert(binPointers[i] == 0);
#endif
    for (unsigned int i=TD->startVertex; i<TD->endVertex; i++)
    {
        prVal = G->attr[i]/outDeg[i];
        for (unsigned int j=G->VI[i]; j<G->VI[i+1]; j++)
        {
            destId = G->EI[j];
            destBin = (destId >> binOffsetBits);
            loc = binPointers[destBin];
            updateBins[destBin][loc] = prVal;
            destIdBins[destBin][loc] = destId;
            binPointers[destBin] = loc+1; 
        }
    }
    for (unsigned int i=0; i<NUM_BINS; i++)
        binPointers[i] = 0;
}

void scatterDPB(threadData* TD, float** updateBins, float** updateBuffers, unsigned int* binPointers, unsigned int* bufferPointers, unsigned int* outDeg)
{
    graph* G = TD->G;
    float prVal = 0.0;
    unsigned int destId = 0;
    unsigned int destBin = 0;
    unsigned int loc = 0;
#ifdef DEBUGL2
    assert((bufferSize%streamSize)==0);
    for (unsigned int i=0; i<NUM_BINS; i++)
    {
        assert(binPointers[i] == 0);
        assert(bufferPointers[i] == 0);
    }
#endif
    for (unsigned int i=TD->startVertex; i<TD->endVertex; i++)
    {
        prVal = G->attr[i]/outDeg[i];
        for (unsigned int j=G->VI[i]; j<G->VI[i+1]; j++)
        {
            destBin = (G->EI[j] >> binOffsetBits);
            updateBuffers[destBin][bufferPointers[destBin]++] = prVal;
            if (bufferPointers[destBin]==bufferSize)             
            {
//                __assume_aligned(updateBuffers[destBin], 64);
                for (unsigned int k=0; k<bufferSize; k+=8)
                {
                    __m256 dummy = _mm256_setr_ps(updateBuffers[destBin][k],updateBuffers[destBin][k+1],updateBuffers[destBin][k+2],updateBuffers[destBin][k+3],updateBuffers[destBin][k+4],updateBuffers[destBin][k+5],updateBuffers[destBin][k+6],updateBuffers[destBin][k+7]);
                    _mm256_stream_ps((updateBins[destBin]+binPointers[destBin]), dummy);
                    binPointers[destBin] += streamSize;
                }
                bufferPointers[destBin] = 0;
            }
        }
    }
    for (unsigned int i=0; i<NUM_BINS; i++)
    {
        for (unsigned int j=0; j<bufferPointers[i]; j++)
            updateBins[i][binPointers[i]++] = updateBuffers[i][j];
    }
    for (unsigned int i=0; i<NUM_BINS; i++)
    {
        binPointers[i] = 0;
        bufferPointers[i] = 0;
    }
}


//////////////////////////////////////////
//// for separate bins to each thread ////
//////////////////////////////////////////

///*****************************************
void gather(graph* G, float*** updateBins, unsigned int*** destIdBins, unsigned int binId, unsigned int** binSizes)
{
    unsigned int startVertex = (binId << binOffsetBits);
    unsigned int endVertex = ((binId+1) << binOffsetBits);
    endVertex = (endVertex > G->numVertex) ? G->numVertex : endVertex;
    for (unsigned int i=startVertex; i<endVertex; i++)
        G->attr[i] = 0;
    
    double tempVal = 0.0;
   
    for (unsigned int i=0; i<NUM_THREADS; i++)
    {
        for (unsigned int j=0; j<binSizes[i+1][binId]; j=j+1) 
        {
#ifdef DEBUGL2
        assert(destIdBins[i][binId][j] < endVertex);
        assert(destIdBins[i][binId][j] >= startVertex);
#endif
        

        G->attr[destIdBins[i][binId][j]] += updateBins[i][binId][j];

        }
    }

    for (unsigned int i=startVertex; i<endVertex; i++)
        G->attr[i] = dampingFactor + (1-dampingFactor)*G->attr[i];
}
//*****************************************/
