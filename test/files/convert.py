f = open('samples_id_tb.txt', 'r')
w = open('samples_fastas_tb.txt', 'w')

for line in f:
	idstr = str(line).strip()
	lineo = "{0}\t/mnt/microbio/ndm-hicf/ogre/pipeline_output/{0}/MAPPING/2e6b7bc7-f52c-4649-8538-c984ab3894bb_R00000039/STD/basecalls/{0}_v3.fasta.gz\n".format(idstr)
	w.write(lineo)
