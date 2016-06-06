#!/bin/bash

echo "" > statistics_saur.txt

for fastas in 10 40 100 400 1000 4000
do
	for threads in 16 8 4 2 1
	do
		echo "threads=$threads, fastas=$fastas" >> statistics_saur.txt
		var=$( { /usr/bin/time -f ", %e real sec, %M mem kb" ./bugmat -t $threads -s test/files/samples_fastas_saur_$fastas.txt -o test/fastas_saur_$fastas; } 2>&1 )
		echo $var >> statistics_saur.txt
	done
done
