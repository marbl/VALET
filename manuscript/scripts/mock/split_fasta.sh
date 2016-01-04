#!/usr/bin/bash
## Script to split a fasta file with multiple sequences into a single file

while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
        outfile=${line:1:20}.fa
        echo $line > $outfile
    else
        echo $line >> $outfile
    fi
done < reference.refseq.fna
