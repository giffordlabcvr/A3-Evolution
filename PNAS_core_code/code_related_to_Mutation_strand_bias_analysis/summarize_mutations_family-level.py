#!/usr/bin/env python

"""
#Description
Summarize mutations in each TE family

#Usage
python summarize_mutations_family-level.py \
       <mutation info in each TE integrant> \
       > <output file>

#Output
1, class: class
2, family: family
3, a.site.q: the number of A sites in the consensus sequence (total value among all alignments)
4, t.site.q: the number of T sites in the consensus sequence (total value among all alignments)
5, g.site.q: the number of G sites in the consensus sequence (total value among all alignments)
6, c.site.q: the number of C sites in the consensus sequence (total value among all alignments)
7, ac: the number of A-to-C mutations (total value among all alignments)
8, ag: the number of A-to-G mutations (total value among all alignments)
9, at: the number of A-to-T mutations (total value among all alignments)
10, ca: the number of C-to-A mutations (total value among all alignments)
11, cg: the number of C-to-G mutations (total value among all alignments)
12, ct: the number of C-to-T mutations (total value among all alignments)
13, ga: the number of G-to-A mutations (total value among all alignments)
14, gc: the number of G-to-C mutations (total value among all alignments)
15, gt: the number of G-to-T mutations (total value among all alignments)
16, ta: the number of T-to-A mutations (total value among all alignments)
17, tc: the number of T-to-C mutations (total value among all alignments)
18, tg: the number of T-to-G mutations (total value among all alignments)
"""


import sys
import numpy as np
argvs = sys.argv

f = open(argvs[1])
sw_thresh = 1000
header_l = f.next().strip().split()
header_l = header_l[6:]

count_sum_d = {}
for line in f:
  line = line.strip().split()
  TE_locus = line[0]
  sw_score = int(line[4])
  if sw_score > sw_thresh:
    nt_compos_l = line[6:10]
    mut_freq_l = line[10:]
    nt_compos_l = [int(c) for c in nt_compos_l]
    mut_freq_l = [int(c) for c in mut_freq_l]
    mut_freq_a = np.array(mut_freq_l)
    nt_compos_a = np.array(nt_compos_l)
    Class,family = TE_locus.split("|")[:2]
    type_t = (Class,family)
    if type_t not in count_sum_d:
      count_sum_d[type_t] = [np.array([0]*4),np.array([0]*12)]
    count_sum_d[type_t][0] = count_sum_d[type_t][0] + nt_compos_a
    count_sum_d[type_t][1] = count_sum_d[type_t][1] + mut_freq_a


print "class\tfamily\t" + "\t".join(header_l)
for type_t in count_sum_d:
  Class,family = type_t
  sum_nt_compos_l = list(count_sum_d[type_t][0])
  sum_mut_freq_l = list(count_sum_d[type_t][1])
  res_l = [Class,family] + sum_nt_compos_l + sum_mut_freq_l
  print "\t".join([str(c) for c in res_l])

