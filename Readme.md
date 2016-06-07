# Abstract
Using nucleotide sequence information, stored in a series of FASTA files, BugMat generates a pairwise distance matrix.  

## Requirements
BugMat is a C++ command line executable.

It can be compiled on Linux (using gcc) and on Windows (using MinGW).

It has few dependencies; it uses the C++ Standard Library 14.
OpenMP is required for parallisation.  

## Performance
BugMat builds in-memory distance matrices.  Using 16 cores on a machine equipped with Intel Xeon E5-2680-v2 Processors (2.8GHz):
* Building a matrix from 400 _M. tuberculosis_ genomes (4.4 million bases per sequence) required about 2GB RAM and took 100 seconds.
* Building a matrix from 4,000 _M. tuberculosis_ genomes (4.4 million bases per sequence) required 20GB RAM, and took about 12,000 seconds.

### Prepare the computer
First of all you should check if the system has gcc compiler and openmp library, the examples are in linux not windows
1- Check openMP library: echo | cpp -fopenmp -dM | grep -i open
2- Check gcc compiler: gcc --version ## get compiler version

Install gcc:
sudo apt-get install gcc-4.2

Install openmp:
apt-get install libgomp1

http://openmp.org/wp/openmp-compilers/
https://huseyincakir.wordpress.com/2009/11/05/installing-openmp-in-linux-debian/

## To compile:
>make clean
>make

or

>g++ -std=c++11 -fopenmp -O3 bugmat.cpp -lz -o bugmat

## Running BugMat
  > ./bugmat --threads 8 --samples file.txt --output folder

### Options:
*  --h, --help                          show this help message and exit
*  --t, --threads numthreads            number of threads for openmp library
*  --s, --sample samplefoldername       path to set of samples sequences
*  --o, --output folder                 folder name for the output files

### Example:
./bugmat -t 8 -s fastas/file_fastas.txt -o fastas
