convert edge list text file to csr binary file

Steps to run:
1. make
2. ./a.out \<inputGraphFileName\> \<outputFileName\>

Input file should be a list of (src,dst) pairs (see exGraph.txt)

For generating weighted graphs,

1. if the input file has weights, use <br />
    ./a.out \<inputGraphFileName\> \<outputFileName\> -w 0 <br />
    **Input file should be a list of (src,dst,wt) tuples (see exGraphWt.txt)

2. if all weights are 1, use <br />
    ./a.out \<inputGraphFileName\> \<outputFileName\> -w 1 <br />
    **Input file should be a list of (src,dst) pairs (see exGraph.txt)

3. if weights are to be randomly assigned (between 1 and 10), use <br />
    ./a.out \<inputGraphFileName\> \<outputFileName\> -w 2 <br />
    **Input file should be a list of (src,dst) pairs (see exGraph.txt)


To create CSC or transpose graph (for pull direction computation), use
1. ./a.out \<inputGraphFileName\> \<outputFileName\> -r \<transposeFileName\>
