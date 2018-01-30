#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#include "../../include/graph.hpp"
#ifndef SORT_HEADER
#include "../../include/sort.hpp"
#endif

using namespace std;

typedef struct threadData
{
    int tid;
    int startVertex;
    int endVertex;
    graph* G;
    unsigned int* outDeg;
    float* srcArr;
    float* dstArr;
}threadData;

#define DEBUG
#undef DEBUG

#define DUMP
//#undef DUMP

#define PERF_MON
#undef PERF_MON

unsigned int MAX_ITER=20;


int NUM_THREADS = 16;
float d = 0.15;
pthread_barrier_t barrier;

//function for openmp
void prOmp (threadData* TD);


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
    
    int rc;
    // graph object
    graph G;

    // read csr file
    if (read_csr(argv[1], &G)==-1)
    {
        printf("couldn't read %s\n", argv[1]);
        exit(1);
    }

#ifdef DEBUG
    printf("%d %d\n", G.numVertex, G.numEdges);
#endif    
    // output Degree array
    float* temp = new float[G.numVertex]();
    unsigned int* outDeg = new unsigned int [G.numVertex]();
    findOutDeg(&G, outDeg);
#ifdef DEBUG
    printf("out degree computed\n");
#endif



    // initialize page rank attribute to 1
    initGraphPullPR (&G, outDeg);
#ifdef DEBUG
    printf("page rank values initialized\n");
#endif 

    //compute transpose for pull based PR
    transposeCSR (&G);
#ifdef DEBUG
    printf("graph transpose computed for pull baseline\n");
#endif 

    //thread attributes and objects
    threadData* TD = (threadData*) malloc (sizeof(threadData)*NUM_THREADS);

    unsigned int numEdgesPerThread = G.numEdges/NUM_THREADS;

    int vcount = -1;
    for (unsigned int i=0; i<NUM_THREADS; i++)
    {
        TD[i].tid = i;
        TD[i].G = &G;
        TD[i].outDeg = outDeg;
        TD[i].startVertex = vcount + 1;
        TD[i].srcArr = temp;
        TD[i].dstArr = G.attr;
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


	struct timespec start, end; 
	float time;
	// measure the start time here
    unsigned int numIter = 0;
    ///////////////////////////////////
    ///// warmup iteration/////////////
    //////////////////////////////////
	if( clock_gettime(CLOCK_REALTIME, &start) == -1) { perror("clock gettime");}
    #pragma omp parallel for schedule(static, 1) num_threads(NUM_THREADS)
    for (unsigned int i=0; i<NUM_THREADS; i++)
    {
        if (numIter%2)
        {
            TD[i].srcArr = temp;
            TD[i].dstArr = G.attr;
            prOmp(&TD[i]);
        }
        else
        {
            TD[i].dstArr = temp;
            TD[i].srcArr = G.attr;
            prOmp(&TD[i]);
        }
    }
    numIter++;

    while(numIter < MAX_ITER+1)
    {
        #pragma omp parallel for schedule(static, 1) num_threads(NUM_THREADS)
        for (unsigned int i=0; i<NUM_THREADS; i++)
        {
            if (numIter%2)
            {
                TD[i].srcArr = temp;
                TD[i].dstArr = G.attr;
                prOmp(&TD[i]);
            }
            else
            {
                TD[i].dstArr = temp;
                TD[i].srcArr = G.attr;
                prOmp(&TD[i]);
            }
        }
        numIter++;
    }


    if (numIter%2==0)
    {
        for (unsigned int i=0; i<G.numVertex; i++)
        {
            if (outDeg[i]>0)
                G.attr[i] = G.attr[i]*outDeg[i];
        }
    }
    else
    {
        for (unsigned int i=0; i<G.numVertex; i++)
        {
            if (outDeg[i]>0)
                temp[i] = temp[i]*outDeg[i];
        }
    }


	if( clock_gettime( CLOCK_REALTIME, &end) == -1 ) { perror("clock gettime");}		
	time = (end.tv_sec - start.tv_sec)+ (int)(end.tv_nsec - start.tv_nsec)/1e9;
    printf("%s, time = %lf\n ", argv[1], time);



#ifdef DUMP
    FILE* fdump = fopen("dumpPR.txt", "w");
    if (fdump == NULL)
    {
        fputs("file error\n", stderr);
        exit(1);
    }
	unsigned int printVertices = (G.numVertex>1000) ? 1000 : G.numVertex;
	if (numIter%2==0)
	{
	    mergeSortWOkey<float>(G.attr, 0, G.numVertex-1);
		for (unsigned int i=0; i<printVertices; i++)
			fprintf(fdump, "%lf\n", G.attr[i]);
	}
	else
	{
	    mergeSortWOkey<float>(temp, 0, G.numVertex-1);
		for (unsigned int i=0; i<printVertices; i++)
			fprintf(fdump, "%lf\n", temp[i]);
	}
    fclose(fdump);
#endif
    //////////////////////////////////////

    // free allocated memory//
    freeMem(&G);
    free(TD);
    delete[] temp;
    delete[] outDeg;
    

    return 0;
}


void prOmp (threadData* TD)
{
    graph* G = TD->G;
    unsigned int* degArr = TD->outDeg;
    float* srcArr = TD->srcArr;
    float* dstArr = TD->dstArr;
    unsigned int start = TD->startVertex;
    unsigned int end = TD->endVertex;
    end = (end > G->numVertex) ? G->numVertex : end;
    float tempVal = 0; 
    unsigned int iter=0;
    for (unsigned int i=start; i<end; i++)
    {
        float tempVal = 0; 
        for (unsigned int j=G->VI[i]; j<G->VI[i+1]; j++)
            tempVal += (srcArr[G->EI[j]]);
		dstArr[i] = (d + (1-d)*tempVal);
        if (degArr[i] > 0)
            dstArr[i] = dstArr[i]/degArr[i];
    }
}

