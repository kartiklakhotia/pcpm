#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#define SORT_HEADER

#define LSB_MASK 0x00000001

unsigned int computeMean (unsigned int a, unsigned int b)
{
    unsigned int meanVal = (a >> 1) + (b >> 1) + ((a & LSB_MASK) + (b & LSB_MASK))/2;
    return meanVal;
}

template <typename T1, typename T2>
void mergeAscending (T1* arr, T2* key, unsigned int low, unsigned int high, unsigned int mid)
{
    unsigned int i = low;
    unsigned int j = mid+1;
    unsigned int k = 0;
    T1* temp = new T1 [high-low+1];
    T2* tempKey = new T2 [high-low+1];
    while(i<=mid && j<=high)
    {
        if (arr[i] <= arr[j])
        {
            tempKey[k] = key[i];
            temp[k++] = arr[i++];
        }
        else
        {
            tempKey[k] = key[j];
            temp[k++] = arr[j++];
        }
    }
    while(i<=mid)
    {
        tempKey[k] = key[i];
        temp[k++] = arr[i++];
    }
    while(j<=high)
    {
        tempKey[k] = key[j];
        temp[k++] = arr[j++];
    }
    for (i=low; i<=high; i++)
    {   
        key[i] = tempKey[i-low];
        arr[i] = temp[i-low];
    }
    delete[] temp;
    delete[] tempKey;
}

template <typename T1, typename T2>
void mergeSortAscending (T1* arr, T2* key, unsigned int low, unsigned int high)
{

    if (low >= high)
        return;
    if ((high - low) == 1)
    {
        if (arr[high] < arr[low])
        {
            T1 temp = arr[high];
            arr[high] = arr[low];
            arr[low] = temp;
            T2 temp2 = key[high];
            key[high] = key[low];
            key[low] = temp2;
        }
        return;
    }
    unsigned int mid = computeMean(low, high);
    mergeSortAscending<T1, T2>(arr, key, low, mid);
    mergeSortAscending<T1, T2>(arr, key, mid + 1, high);
    mergeAscending<T1, T2>(arr, key, low, high, mid);
    
} 

template <typename T1, typename T2>
void merge (T1* arr, T2* key, unsigned int low, unsigned int high, unsigned int mid)
{
    unsigned int i = low;
    unsigned int j = mid+1;
    unsigned int k = 0;
    T1* temp = new T1 [high-low+1];
    T2* tempKey = new T2 [high-low+1];
    while(i<=mid && j<=high)
    {
        if (arr[i] >= arr[j])
        {
            tempKey[k] = key[i];
            temp[k++] = arr[i++];
        }
        else
        {
            tempKey[k] = key[j];
            temp[k++] = arr[j++];
        }
    }
    while(i<=mid)
    {
        tempKey[k] = key[i];
        temp[k++] = arr[i++];
    }
    while(j<=high)
    {
        tempKey[k] = key[j];
        temp[k++] = arr[j++];
    }
    for (i=low; i<=high; i++)
    {   
        key[i] = tempKey[i-low];
        arr[i] = temp[i-low];
    }
    delete[] temp;
    delete[] tempKey;
}

template <typename T1, typename T2>
void mergeSort (T1* arr, T2* key, unsigned int low, unsigned int high)
{

    if (low >= high)
        return;
    if ((high - low) == 1)
    {
        if (arr[high] > arr[low])
        {
            T1 temp = arr[high];
            arr[high] = arr[low];
            arr[low] = temp;
            T2 temp2 = key[high];
            key[high] = key[low];
            key[low] = temp2;
        }
        return;
    }
    unsigned int mid = computeMean(low, high);
    mergeSort<T1, T2>(arr, key, low, mid);
    mergeSort<T1, T2>(arr, key, mid + 1, high);
    merge<T1, T2>(arr, key, low, high, mid);
    
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

