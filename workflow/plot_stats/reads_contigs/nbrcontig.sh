#!/bin/bash
touch 'contigs.csv'
for file in *_reads_length.json
do
	basename $file _reads_length.json >> 'contigs.csv'
	jq -r 'to_entries | .[-1].key' $file >> 'contigs.csv'
done
