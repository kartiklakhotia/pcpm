convert edge list text file to csr binary file

Steps to run:
1. make
2. ./a.out \<inputGraphFileName\> \<outputFileName\>


Input file should have a list of (src,dst) pairs.

To create CSC (for pull direction computation), use createReverseCSR function
in graph.h
