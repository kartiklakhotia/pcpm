convert edge list text file to csr binary file

Steps to run:
1. make
2. ./a.out \<inputGraphFileName\> \<outputFileName\>

Input file should be a list of (src,dst) pairs (see exGraph.txt)

For weighted graphs, use
./a.out \<inputGraphFileName\> \<outputFileName\> -w
Input file should be a list of (src,dst) pairs (see exGraphWt.txt)


To create CSC (for pull direction computation), use createReverseCSR function
in graph.h
