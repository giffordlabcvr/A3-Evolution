#!/usr/bin/env bash

#select seq
python select_alignment_site.py \
  demo_data/original_MSA.fas \
  demo_data/second_MSA \
  0.85 0

#extract branch_info
python remove_outlierOTUs.py \
       demo_data/tree_first.nwk \
       demo_data/second_MSA_seq_selected_pos.fas \
       demo_data/rm_outliers_tree_first \
       5

