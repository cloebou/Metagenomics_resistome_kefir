#!/bin/bash
touch 'read_classified_kraken.csv'
for file in */*.kraken 
do 
	awk '{print $2}' $file | awk -F '.' 'END{print}' >> 'read_classified_kraken.csv'
done
