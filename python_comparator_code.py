## Script that takes in a list of fastas and generates  summaries of their differences
# this is python implementation to which bugMat was compared in our publication.

# code by Drs Madeleine Cule and David Eyre, Oxford University

import sys, re, os, gzip
from Bio import SeqIO
from itertools import izip

bases = "ACGT"

def get_distance( seq1, seq2 ):
	'''Function to calculate the [Hamming] distance between two sequences'''
	return sum(True for c1, c2 in izip( seq1, seq2 ) if c1!=c2 and c1 in bases and c2 in bases )


if __name__=="__main__":

	## Parse command line arguments
	listoffasta, outname_prefix = sys.argv[1:]  
	
	outname = outname_prefix
	## Read in all the sequences, and replace the id with the required nicename
	seqlist = []
	with open( listoffasta ) as fp:
		for line in fp:
			nicename, fapath = line.strip().split()
			if (os.path.exists( fapath)):
				f = gzip.open( fapath )
				fa = SeqIO.read( f, "fasta" )
				fa.id = nicename
				seqlist.append( fa )
			else:
				sys.stderr.write(fapath+" does not exist, skipping...")
	
	## Find nonshared positions
	seq_generator = izip( *seqlist )
	nonshared_pos =[ i for ( i, a ) in enumerate( seq_generator ) 
					   if len( set( [ ai for ai in a if ai in bases ] ) ) >  1 ]
	sys.stderr.write("Successfully obtained nonshared_diffs; there are %s of them.\n"%len( nonshared_pos ) )
	
	## Sort out the SNPs
	for seq in seqlist:
		nonshared_bases = "".join( seq.seq[ i ] for i in nonshared_pos )
		seq.seq._data = nonshared_bases
	SeqIO.write( seqlist , "%s_snps.fa"%outname, "fasta" )
	sys.stderr.write("Successfully wrote snps fasta file.\n")
	
	## Write the positions
	with open("%s_positions.txt"%outname, "w" ) as out:
		out.write( "\n".join( [ str( n+1 ) for n in nonshared_pos ] ) )
		out.write("\n")
	sys.stderr.write("Successfully wrote nonshared positions.\n")
	
	## Do the data matrix
	mfasta = "%s_snps.fa"%outname
	listofsamples = []
	
	strings = dict()
	for seq_record in SeqIO.parse( mfasta, "fasta" ):
		strings[ seq_record.id ] = str( seq_record.seq )
		listofsamples.append( seq_record.id )
	
	listofsamples2 = list(listofsamples)
	
	with open("%s.dat"%outname, "w" ) as fp:
		fp.write("%s\n"%("\t".join( listofsamples ) ) ),
		for s1 in listofsamples:
			fp.write( "%s\t"%s1 )
			for s2 in listofsamples2:
				fp.write( "%s\t"%get_distance( strings[ s1 ], strings[ s2 ] ) )
			fp.write("\n")
		fp.close()
	
	sys.stderr.write("Successfully wrote pairwise distance matrix.\n")
	sys.stderr.write("Done.")
