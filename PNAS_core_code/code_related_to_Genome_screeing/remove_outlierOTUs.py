#!/usr/bin/env python

"""
#Description
Remove OTUs having an extraordinary long external branch

#Usage
python ./remove_outlierOTUs.py \
       <first tree (nwk format)> \
       <original multiple sequence alignment file (multi-fasta format)> \
       <output file prefix> \
       <max Z score of the external branch length>

#Output
1, <output file prefix>_bl.info.txt
Information of 1) OTU name, 2) external branch length, and 3) Z score 

2, <output file prefix>_rm_outlier.fas
Multile sequence alignment (multi-fasta format) without the sequences of outlier-OTUs

3, <output file prefix>_outlier.fas
Multile sequence alignment for the sequences of outlier-OTUs
"""


import sys, re
import numpy as np
argvs = sys.argv

tree_f = open(argvs[1])
seq_f = open(argvs[2])
out_prefix = argvs[3]
thresh_z_score = float(argvs[4]) 

tree_line = tree_f.next().strip()

bl_d = {}
for m in re.finditer("[0-9]+:[0-9\.]+", tree_line):
  otu_bl = m.group()
  otu = re.sub(r'([0-9]+):.+',r'\1',otu_bl)
  bl = re.sub(r'[0-9]+:(.+)',r'\1',otu_bl)
  bl = float(bl)
  otu = otu.replace(' ','')
  bl_d[otu] = bl

bl_a =  np.array(bl_d.values())
ave = np.average(bl_a)
sd = np.std(bl_a)
thresh = ave + sd * thresh_z_score

seq_f = open(argvs[2])
seq_d = {}
for line in seq_f:
  line = line.strip()
  if '>' in line:
    otu = line.replace('>','')
    otu = otu.replace('_','')
    seq_d[otu] = ""
  else:
    seq = line
    seq_d[otu] = seq_d[otu] + seq

#bl_info
out_bl_info_f_name = out_prefix + '_bl.info.txt'

out_bl_info_f = open(out_bl_info_f_name,"w")
print >> out_bl_info_f, 'OTU\texternal_branch_len\tz_score'
for otu in bl_d:
  bl = bl_d[otu]
  z_score = (bl - ave) / sd
  res_l = [otu,bl,z_score]
  print >> out_bl_info_f, "\t".join([str(c) for c in res_l]) 


#fasta

out_fasta_f_name = out_prefix + '_rm_outlier.fas'
out_fasta__excluded_f_name = out_prefix + '_outlier.fas'

out_fasta_f = open(out_fasta_f_name,"w")
out_fasta__excluded_f = open(out_fasta__excluded_f_name,"w")

for otu in seq_d:
  bl = bl_d[otu]
  seq = seq_d[otu]
  res = '>' + otu + "\n" + seq
  if bl < thresh:
    print >> out_fasta_f, res
  else:
    print >> out_fasta__excluded_f, res

