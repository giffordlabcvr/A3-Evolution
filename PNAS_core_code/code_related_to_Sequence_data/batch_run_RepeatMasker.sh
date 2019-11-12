#!/usr/bin/env bash
#RepeatMasker version open-4.0.9

################Usage################
# bash batch_run_RepeatMasker.sh \
#          <species name in NCBI Taxonomy Database> \
#          <whole genome sequence fasta file>

#args
species=$1
genome_file=$2

RepeatMakser -q -xsmall -a -pa 4 -species ${species} ${genome_file} 

