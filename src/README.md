### Directory structure ###

* pdpr -> pull direction PageRank
* propagationBlocking-> vertex centric gas with streaming stores (bypasses cache and reduces communication)
* pcpm -> partition centric gas

### Steps to compile and run ###
1. cd <directory>
2. make
3. ./pr <csr_binary_dataset> <NUM_THREADS> <NUM_ITERATIONS>
