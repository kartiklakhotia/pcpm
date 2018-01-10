#include "graph.h"

///SORT functions to sort edges on source or destination//
///These functions are internal. Not to be exposed outside
#define LSB_MASK 0x00000001

unsigned int computeMean (unsigned int a, unsigned int b)
{
    unsigned int meanVal = (a >> 1) + (b >> 1) + ((a & LSB_MASK) + (b & LSB_MASK))/2;
    return meanVal;
}

template <typename T>
void mergeWOkey (T* arr, unsigned int low, unsigned int high, unsigned int mid)
{
    unsigned int i = low;
    unsigned int j = mid+1;
    unsigned int k = 0;
    T* temp = new T [high-low+1];
    while(i<=mid && j<=high)
    {
        if (arr[i] >= arr[j])
            temp[k++] = arr[i++];
        else
            temp[k++] = arr[j++];
    }
    while(i<=mid)
        temp[k++] = arr[i++];
    while(j<=high)
        temp[k++] = arr[j++];
    for (i=low; i<=high; i++)
        arr[i] = temp[i-low];
    delete[] temp;
}


template <typename T>
void mergeSortWOkey (T* arr, unsigned int low, unsigned int high)
{

    if (low >= high)
        return;
    if ((high - low) == 1)
    {
        if (arr[high] > arr[low])
        {
            T temp = arr[high];
            arr[high] = arr[low];
            arr[low] = temp;
        }
        return;
    }
    unsigned int mid = computeMean(low, high);
    mergeSortWOkey<T>(arr, low, mid);
    mergeSortWOkey<T>(arr, mid + 1, high);
    mergeWOkey<T>(arr, low, high, mid);
    
} 

// End of sort functions //
///////////////////////////

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


    G->VI = new unsigned int[G->numVertex];

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

    G->EI = new unsigned int[G->numEdges];
    fread (G->EI, sizeof(unsigned int), G->numEdges, graphFile);
    if (feof(graphFile))
    {
        delete[] G->VI;
        delete[] G->EI;
        printf("unexpected end of file while reading edges\n");
        return -1;
    }
    else if (ferror(graphFile))
    {
        delete[] G->VI;
        delete[] G->EI;
        printf("error reading file\n");
        return -1;
    }
    fclose(graphFile);
    return 1;
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

void printGraph (graph* G)
{
    printf("number of vertices are %d\n", G->numVertex);
    printf("number of edges are %d\n", G->numEdges);
    for (unsigned int i=0; i<G->numVertex; i++)
        printf("%d ", G->VI[i]);
    printf("\n");
    for (unsigned int i=0; i<G->numEdges; i++)
        printf("%d ", G->EI[i]);
    printf("\n");
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

void createReverseCSR(graph* G1, graph* G2)
{
    G2->numVertex = G1->numVertex;
    G2->numEdges = G1->numEdges;
    G2->VI = new unsigned int[G1->numVertex]();
    G2->EI = new unsigned int[G1->numEdges]; 

    for (unsigned int i=0; i<G1->numEdges; i++)
    {
        if (G1->EI[i] < G1->numVertex-1)
            G2->VI[G1->EI[i]+1]++;
    }
    for (unsigned int i=1; i<G1->numVertex; i++)
        G2->VI[i] += G2->VI[i-1];
    unsigned int* tempId = new unsigned int [G1->numVertex]();
    for (unsigned int i=0; i<G1->numVertex; i++)
    {
        unsigned int maxId = (i==G1->numVertex-1) ? G1->numEdges : G1->VI[i+1];
        for (unsigned int j=G1->VI[i]; j<maxId; j++)
        {
            G2->EI[G2->VI[G1->EI[j]] + tempId[G1->EI[j]]] = i;
            tempId[G1->EI[j]]++;
        } 
    }
    delete[] tempId;
}


void freeMem (graph* G)
{
    delete[] G->VI;
    delete[] G->EI;
}



