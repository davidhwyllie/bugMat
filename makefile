all: program

program: 
	g++ -std=c++11 -fopenmp -O3 bugmat.cpp -lz -o bugmat

clean:
	rm -rf *o bugmat
