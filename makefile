all: plain

compressed: 
	g++ -std=c++14 -fopenmp -O3 bugmat.cpp -lz -o bugmat_compressed

plain:
	g++ -std=c++14 -fopenmp -O3 bugmat_on_fasta.cpp -o bugmat_plain

clean:
	rm -rf *o bugmat_compressed bugmat_plain
