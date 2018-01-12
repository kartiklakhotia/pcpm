Greedy labeling code. To run, enter the following commands:

* make
* ./relabel <inputFileName> <outputFileName> <traversalOrder> <dumpLabels>



Following arguments need to be specified:

* inputFileName
* outputFileName



Following arguments are optional:

* traversalOrder - the order in which nodes are traversed. Default value = 3.
                   (0 -> original labels; 1 -> increasing order of degree; \
                    2 -> decreasing order of degree; 3 -> randomized)

* dumpLabels - output a text file with new labels for each node. Default value = 0.
               (no output if value is 0; output relabel.txt otherwise)
