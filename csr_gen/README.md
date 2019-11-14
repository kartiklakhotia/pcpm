### convert edge list text file to csr binary file ###

Steps to run:
1. make
2. ./a.out \<inputGraphFileName\> \<outputFileName\>

Input file should be a list of directed (src,dst) pairs (see exGraph.txt) <br />

If the input file represents an undirected graph, specify using the "-u" option eg. <br />

 ./a.out \<inputGraphFileName\> \<outputFileName\> -u  <br />



### For generating weighted graphs, ###

1. if the input file has weights, use <br />
    ./a.out \<inputGraphFileName\> \<outputFileName\> -w 0 <br />
    **Input file should be a list of (src,dst,wt) tuples (see exGraphWt.txt)

2. if all weights are 1, use <br />
    ./a.out \<inputGraphFileName\> \<outputFileName\> -w 1 <br />
    **Input file should be a list of (src,dst) pairs (see exGraph.txt)

3. if weights are to be randomly assigned (between 1 and 10), use <br />
    ./a.out \<inputGraphFileName\> \<outputFileName\> -w 2 <br />
    **Input file should be a list of (src,dst) pairs (see exGraph.txt)


### Converting Large Graphs ### 
Enable the following flags in makefile(s):
1. HUGE\_EDGE: if number of edges is >4B
2. HUGE\_VERTEX: if number of vertices is >2B 


### To create CSC or transpose graph (for pull direction computation) ###
1. ./a.out \<inputGraphFileName\> \<outputFileName\> -r \<transposeFileName\>
