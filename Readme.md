# Abstract
Using nucleotide sequence information, stored in a series of FASTA files, BugMat generates a pairwise distance matrix.  

## Requirements
BugMat is a C++ command line executable.

It can be compiled on Linux and on Windows.

It uses the C++ Standard Library 14 and the zlib library.
OpenMP is required for parallelisation.  

## Performance
BugMat builds in-memory distance matrices.  Using 16 cores on a machine equipped with Intel Xeon E5-2680-v2 Processors (2.8GHz):
* Building a matrix from 400 _M. tuberculosis_ genomes (4.4 million bases per sequence) required about 2GB RAM and took 100 seconds.
* Building a matrix from 4,000 _M. tuberculosis_ genomes (4.4 million bases per sequence) required 20GB RAM, and took about 12,000 seconds.

### Prepare the computer
First of all you should check if the system has gcc compiler and openmp library, the examples are in linux not windows:  

1- Check openMP library: 
```
echo | cpp -fopenmp -dM | grep -i open
```

2- Check gcc compiler:   
```
gcc --version ## get compiler version
```

Install gcc:
```
sudo apt-get install gcc-4.2
```

Install openmp:
```
apt-get install libgomp1
```

cf. 
http://openmp.org/wp/openmp-compilers/
https://huseyincakir.wordpress.com/2009/11/05/installing-openmp-in-linux-debian/

## To compile:
>make clean
>make

or

>g++ -std=c++11 -fopenmp -O3 bugmat.cpp -lz -o bugmat

## Running BugMat
  > ./bugmat --threads 8 --samples file.txt --output folder

### Compilation on Windows
We have successfully compiled BugMat on windows using both DevC++ and MS Visual Studio 15.
We have tested this on Windows 10/8/7 systems.
Running the software on Windows 7 is complicated by a a known issue with zlib 1.2.8 on Windows 7 (but not later systems) which is described, and can be addressed as described, here:
http://www.tannerhelland.com/5076/compile-zlib-winapi-wapi-stdcall/
https://sourceforge.net/p/globalplatform/code/HEAD/tree/trunk/zlib-1.2.8/
https://sourceforge.net/projects/globalplatform/files/zLib/


### Options:
*  --h, --help                          show this help message and exit
*  --t, --threads numthreads            number of threads for openmp library
*  --s, --sample samplefoldername       path to set of samples sequences
*  --o, --output folder                 folder name for the output files

### Example:
./bugmat -t 8 -s fastas/file_fastas.txt -o fastas


### Expected output
BugMat reports whether each of four stages it undertakes succeeded (0) or failed (1).
If all four stages report zero (success), then the output will be written into the --output folder.

```
# succeeds
> ./bugmat -t 8 -s fastas/file_fastas.txt -o fastas
read&clean,process,write,total
0, 0, 0, 0

# fails, as no output directory  is provided
> ./bugmat -t 8 -s fastas/file_fastas.txt 
read&clean,process,write,total
0, 0, 1, 1
```

### output files
Three output files are produced.

* alignment_samples.fa  
This is a fasta file in which only variant bases are included.

* bugmat_snp.txt
The reports the pairwise SNP distances between all sequences analysed.

* bugmat_count_bases.txt
This is a text representation of the model used by the software to represent the positions of variation.
