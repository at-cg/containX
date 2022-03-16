containX
========================================================================

containX is a prototype implementation of an algorithm that decides which contained reads can be dropped during overlap graph sparsfication. Reads which are substrings of longer reads are typically referred to as contained reads. The string graph model filters out contained reads during graph construction. Contained reads are mostly considered redundant in most assembly algorithms. However, removing all contained reads can lead to coverage gaps, especially in diploid, polyploid genomes and metagenomes. Here we have implemented novel heuristics to distinguish redundant and non-redundant contained reads. 

## Install

Clone source code from master branch. 
  ```sh
  git clone https://github.com/at-cg/containX.git
  ```
To compile, the software requires C++ compiler with c++11 and openmp, which are available by default in GCC >= 4.8.
  ```sh
	cd containX
	make
  ```
Expect `containX` executable in your folder.

## Usage on diploid genomes


