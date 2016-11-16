/**
 * @file bugmat.cpp
 * @name bugmat
 * @author oriol mazariegos
 * @date 29 April 2015
 * @brief Sequences fasta files 1-alignmed data  
 */

#include <time.h>
#include <sys/timeb.h>
#include <omp.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <unordered_map>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <ctime>
#define version "v2016-11-16"

using namespace std;

/************************** FUNCTIONS DECLARATION **************************/
void add_sample(string id, string sample);


/****************************** OPENMP THREADS *****************************/
/** @brief number threads for openmp */
int num_threads = 8;

/********************************* VARIABLES *******************************/
/** @brief max size for each sample */
unsigned long size_sample = 5000000;
/** @brief size of samples */
unsigned long num_samples = 0;
/** @brief list of samples, all in memory */
vector<string> list_samples;
/** @brief variant sites index */
vector<unsigned long> indexes;
/** @brief ref model for invariant sites and N's */
string model;
/** @brief distance matrix result */
vector<vector<unsigned long> > matrix;
/** @brief bases + N */
string bases [5] = {"A","C","T","G","N"};
/** @brief instance structure with the bases, vocabulary of the samples {A,C,G,T,N}  */
vector<vector<unsigned long>> bases_count;

/** @brief ides structures linkage name with index*/
unordered_map<string,unsigned long> ides_name_index;
/** @brief ides structures linkage index with name*/
unordered_map<unsigned long,string> ides_index_name;

/** @brief list of path to samples */
vector<vector<string> > sample_path_list;
/** @brief path to file with list of path to samples */
string input_path = "";
/** @brief output multifasta file */
string output_path = "alignment_samples.fa";
/** @brief output snp file */
string output_snp = "bugmat_snp.txt";
/** @brief output count_bases file */
string output_count_bases = "bugmat_count_bases.txt";
/** @brief output folder */
string output_folder = "";

/**************************** READ/WRITE ***************************/
void read_paths(){
	ifstream myfile (input_path);

	string nice_name;
	string sample_path;
	string line;
	while ( getline (myfile,line) )
	{
		istringstream iss(line);
		iss >> nice_name >> sample_path;
		vector<string> value;
		value.push_back(nice_name);
		value.push_back(sample_path);
		sample_path_list.push_back(value);
	}

	myfile.close();
}


string read_sample(string id,string path){
	char *cstr = new char[path.length() + 1];
	strcpy(cstr, path.c_str());
	
	string sample;
	ifstream myfile;
	myfile.open(cstr);
	myfile >> sample;

	return sample;
}

void write_samples(){
	string path = output_folder+"/"+output_path;
	ofstream myfile (path);

	#pragma omp parallel for num_threads(num_threads)
	for(int i=0; i<list_samples.size(); i++){
		string line;
		string id = ides_index_name[i];
		string sample;
		for(unsigned long k=0; k<indexes.size(); k++){
			unsigned long z = indexes[k];
			sample+=list_samples[i][z];
		}

		line = ">"+id+"\n"+sample+"\n";
		#pragma omp critical
		myfile << line; 
	}

	myfile.close();
}

void write_matrix(){
	string path = output_folder+"/"+output_snp;
	ofstream myfile (path);

	#pragma omp parallel for num_threads(num_threads)
	for(unsigned long i=0; i<num_samples; i++){
		string id1 = ides_index_name[i];
		for(unsigned long j=0; j<num_samples; j++){
			string id2 = ides_index_name[j];
			unsigned long distance = matrix[i][j];	
			#pragma omp critical
			myfile << id1 << "\t" << id2 << "\t" << distance << endl; 
		}	
	}

	myfile.close();
}

void write_count_bases(){
	string path = output_folder+"/"+output_count_bases;
	ofstream myfile (path);

	//print model
	myfile << "model : ";
	for(unsigned long j=0; j<size_sample-1; j++){
		myfile << model[j] << ",";
	}
	myfile << model[size_sample-1] << endl;
	//print bases
	int num_bases = 5;
	int inv_sites;
	for(unsigned long i=0; i<num_bases; i++){
		inv_sites = 0;

		string base_value = bases[i];
		
		myfile << base_value << " : index = ";
		for(unsigned long j=0; j<size_sample-1; j++){
			unsigned long num = bases_count[i][j];
			myfile << num << ",";
		
			if(model.compare(j,1,base_value) == 0){
				inv_sites += 1;	
			}
		}
		unsigned long num = bases_count[i][size_sample-1];
		myfile << num << " : ";

                if(model.compare(size_sample-1,1,base_value) == 0){
                	inv_sites += 1;
                }

		myfile << "inv = " << inv_sites << "";
		myfile << endl;
	}

	myfile.close();
}

/****************************** DISTANCE MATRIX *****************************/
int compare_fast(char a, char b) {

	switch(a) {
		case 'a': a = 'A'; break;
		case 'c': a = 'C'; break;
		case 't': a = 'T'; break;
		case 'g': a = 'G'; break;
	}

	switch(b) {
		case 'a': b = 'A'; break;
		case 'c': b = 'C'; break;
		case 't': b = 'T'; break;
		case 'g': b = 'G'; break;
	}

	bool noerr = 1;
	switch(a) {
		case 'U': case 'u': a = 'T'; break;
		case ' ': case '\t': case '\r': case '\n': case '\v': case '\b': case '\f': noerr = 0; break;
	}
	switch(b) {
		case 'U': case 'u': b = 'T'; break;
		case ' ': case '\t': case '\r': case '\n': case '\v': case '\b': case '\f': noerr = 0; break;
	}
	return (a!=b) && !(a=='N' || b=='N') && !(a=='-' || b=='-') && !(a=='?' || b=='?') && noerr;
}

void process_samples(){
	//initialize matrix
	matrix.resize( num_samples , vector<unsigned long>( num_samples , num_samples ) );

	#pragma omp parallel for num_threads(num_threads)
	for(int i=0; i<num_samples; i++){
		for(int j=i; j<num_samples; j++){
			unsigned long value = 0;
			for(unsigned long k=0; k<indexes.size(); k++){
				unsigned long z = indexes[k];
				value += compare_fast(list_samples[i][z],list_samples[j][z]);
			}
			matrix[i][j] = value;
			matrix[j][i] = value;
		}
	}
}

void clean_data_from_model(){
	indexes.clear();
	for(int i=0; i<size_sample; i++){
		if(model[i] == '.'){
			indexes.push_back(i);
		}
	}
}

/****************************** MODEL *****************************/
void init_model(){
	for(unsigned long j=0; j<size_sample; j++){
		model.append("\t");
	}
}

void change_model(string * sample){
	for(unsigned long i=0; i<size_sample; i++){
		if(model[i] == '\t') {
			model[i] = (*sample)[i];
		}else{
			if(model[i] != '.'){
				if( model[i] == 'N' || model[i] == '-' || model[i] == '?'){
					if(((*sample)[i] != 'N' || (*sample)[i] != '-' || (*sample)[i] != '?')){
						model[i] = (*sample)[i];
					}
				}else if(model[i] != (*sample)[i] && ((*sample)[i]!='N' && (*sample)[i]!='-' && (*sample)[i]!='?')){
					if(!((model[i] == 'U' || model[i] == 'T') && ((*sample)[i] == 'U' || (*sample)[i] == 'T'))){
						model[i] = '.';
					}
				}				
			}
		}
	}
}

/************************* COUNT BASES **************************/
void init_count_bases(){
	vector<unsigned long> A(size_sample);
	vector<unsigned long> C(size_sample);
	vector<unsigned long> T(size_sample);
	vector<unsigned long> G(size_sample);
	vector<unsigned long> N(size_sample);

	bases_count.push_back(A);
	bases_count.push_back(C);
	bases_count.push_back(T);
	bases_count.push_back(G);
	bases_count.push_back(N);

	#pragma omp parallel for num_threads(num_threads)
	for(int i=0; i<size_sample; i++){
		bases_count[0][i]=0;
		bases_count[1][i]=0;
		bases_count[2][i]=0;
		bases_count[3][i]=0;
		bases_count[4][i]=0;
	}
}

void change_count_cases(string * sample){
	#pragma omp parallel for num_threads(num_threads)
	for(int i=0; i<size_sample; i++){
		if((*sample)[i]=='A'){
			bases_count[0][i]++;
		}else if((*sample)[i]=='C'){
			bases_count[1][i]++;
		}else if((*sample)[i]=='T'){
			bases_count[2][i]++;
		}else if((*sample)[i]=='G'){
			bases_count[3][i]++;
		}else{
			bases_count[4][i]++;
		}	
	}
}

/****************************** API *****************************/

void add_sample(string id, string * sample){
	//check the first insertion and fix the length
	if(num_samples == 0){
		size_sample = (*sample).length();
	}

	//add sample to the list
	list_samples.push_back(*sample);

	ides_index_name[num_samples] = id;
	ides_name_index[id] = num_samples;

	num_samples++;
	
	//clean data
	change_model(sample);
	clean_data_from_model();

	//count bases
	change_count_cases(sample);
}

/****************************** PRINTS *****************************/

void print_samples(){
	for(unsigned long i=0; i<list_samples.size(); i++){	
		cout << list_samples[i] << endl;;
	}
}

void print_samples_indexes(){
	for(unsigned long i=0; i<list_samples.size(); i++){	
		for(unsigned long k=0; k<indexes.size(); k++){
			unsigned long z = indexes[k];
			cout << list_samples[i][z];
		}
		cout << endl;
	}
}

void print_indexes(){
	for(unsigned long i=0; i<indexes.size(); i++){	
		cout << indexes[i] << "|";
	}
	cout << endl;
}

void print_model(){
	for(unsigned long i=0; i<size_sample; i++){	
			cout << model[i];
	}
	cout << endl;
}

void print_matrix(){
	for(unsigned long i=0; i<num_samples; i++){
		string id1 = ides_index_name[i];
		for(unsigned long j=0; j<num_samples; j++){
			string id2 = ides_index_name[j];
			unsigned long distance = matrix[i][j];	
			cout << id1 << "\t" << id2 << "\t" << distance << endl; 
		}	
	}
}
/****************************** HELP *****************************/
void print_help(){

	printf("\n");

	printf ("bugmat application HELP - %s\n",version);

	printf("\n");

	printf ("Options:\n");
	printf("  --h, --help,				show this help message and exit\n");
	printf("  --t, --threads numthreads, 		number of threads for openmp library\n");
	printf("  --s, --sample samplefoldername, 	path to set of samples sequences\n");
	printf("  --o, --output folder, 		folder name for the output files\n");

	printf("\n");

	printf("Use starting:\n");
	printf("  > ./name_app --threads 8 --samples file.txt --output folder\n");
	printf("  > ./name_app -t 8 -s file.txt -o folder\n");

	printf("\n");
}
/****************************** INIT *****************************/
void init_system(){
	num_samples = 0;
	init_model();
	init_count_bases();
	clean_data_from_model();
}

void deallocate_memory(){
	sample_path_list.clear();
	list_samples.clear();
	indexes.clear();
	ides_index_name.clear();
}

/****************************** MAIN *****************************/
int main (int argc, char **argv)
{
	int c;
	while (1) {
		static struct option long_options[] =
		{
			//These options don’t set a flag.
			//We distinguish them by their indices. 
			{"help",   	no_argument, 		0, 		'h'},
			{"threads",   	required_argument, 	0, 		't'},
			{"samples",   	required_argument, 	0, 		's'},
			{"output",  	required_argument, 	0, 		'o'},
			{0, 		0, 			0, 		0}
		};
		//getopt_long stores the option index here. 
		int option_index = 0;

		c = getopt_long (argc, argv, "h:t:r:s:o:",long_options, &option_index);

		//Detect the end of the options.
		if (c == -1) { 
			break;
		}

		switch (c) {
			case 'h':
				print_help();
				return 0;
			case 't':
				num_threads = atoi(optarg);
				break;
			case 's':
				input_path = optarg;
				break;
			case 'o':
				output_folder = optarg;
				break;
		}
	}

	init_system();

	string id;
	string path;
	string line;
	//read file with path to fasta files
	read_paths();

	time_t begin;
	time_t end;

	time_t beginloop;
	time_t endloop;

	double total = 0;
	double totaladd = 0;
	double totalread = 0;

	//write head for output time
	cout << "read&clean,process,write,total" << endl;

	//loop for each fasta filea
	time(&beginloop);
	
	#pragma omp parallel for num_threads(num_threads)
	for(int i=0; i<sample_path_list.size(); i++){
		string id = sample_path_list[i][0];
		string path = sample_path_list[i][1];

		//read fasta file
		string sample = read_sample(id, path); 

		//insert sample and create indexes
		#pragma omp critical
		add_sample(id,&sample);
		#pragma omp critical
		sample.clear();
	}
	time(&endloop);
	double totalloop = difftime (endloop,beginloop);
	cout << totalloop << ", ";
	total += totalloop;

	//process distance
	time(&begin);	
	process_samples();	
	time(&end);
        total += difftime (end,begin);
        cout << difftime (end,begin) << ", ";

	//write samples to multifastafile
        time(&begin);
	write_samples();
	//write distance matrix
	write_matrix();
	//write count_bases
	write_count_bases();
        time(&end);
        total += difftime (end,begin);
        cout << difftime (end,begin) << ", ";

        cout << total << endl;

	return 0;
}
