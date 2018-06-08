# README #


### What is this repository for? ###

* Find the paper on https://arxiv.org/abs/1709.07122
* Baseline is propagation blocking (http://www.scottbeamer.net/pubs/beamer-ipdps2017.pdf)
* SpMV acceleration - Avoiding random accesses and reducing overall communication

### Directory Structure ###

* src -> different pagerank implementations (baselines and partition-centric processing)
* include -> hpp files with general purpose functions used in all implementations
 
### Input Data ###

* The programs in src/ directory use graph in csr binary format
* Use the csr_gen utility to create binary csr files from edge list of a graph
