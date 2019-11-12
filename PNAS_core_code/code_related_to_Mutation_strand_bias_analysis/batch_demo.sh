#!/usr/bin/env bash

##Description
#Example to use the codes on the demo data

#calc_mutation_freq
python detect_mutations_in_each_TE_integrant.py \
       demo_data/hg19.fa.head_50000.align \
       > demo_data/mutation_freq_hg19.txt

#sum family
python summarize_mutations_family-level.py \
       demo_data/mutation_freq_hg19.txt \
       > demo_data/mutation_freq_hg19.sum_family.txt

#calc mutation rate
R --vanilla --slave --args \
  demo_data/mutation_freq_hg19.sum_family.txt \
  demo_data/mutation_freq_hg19.sum_family.mut_rate.txt \
  < calc_mut_rate_family-level.R
