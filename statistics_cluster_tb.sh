#!/bin/bash

echo "first execution with cluster samples" > statistics_cluster_tb.txt

for fastas in 10 50 100 200 280
do
	for threads in 16 8 4 2 1
	do
		echo "threads=$threads, fastas=$fastas" >> statistics_cluster_tb.txt
		var=$( { /usr/bin/time -f ", %e real sec, %M mem kb" ./bugmat -t $threads -s test/files/samples_fastas_cluster_tb_$fastas.txt -o test/fastas_cluster_tb_$fastas; } 2>&1 )
		echo  $var >> statistics_cluster_tb.txt
	done
done

echo "second execution with nocluster samples" >> statistics_cluster_tb.txt

for fastas in 10 50 100 200 280
do
	for threads in 16 8 4 2 1
	do
		echo "threads=$threads, fastas=$fastas" >> statistics_cluster_tb.txt
		var=$( { /usr/bin/time -f ", %e real sec, %M mem kb" ./bugmat -t $threads -s test/files/samples_fastas_nocluster_tb_$fastas.txt -o test/fastas_nocluster_tb_$fastas; } 2>&1 )
		echo $var >> statistics_cluster_tb.txt
	done
done
