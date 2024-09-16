#!/usr/bin/env python

"""
#Description
Detect mutations in each TE integrant compared to the consensus sequence

#Usage
python detect_mutations_in_each_TE_integrant.py \
       <RepeatMasker output .algn file> mutation_freq_hg19.txt \
       > <output file>

#Output
1, TE_locus: TE integrant name
2, chrom: chromosome
3, start: start
4, end: end
5, sw_score: SW score
6, total_seq_len: sequence length of the integrant
7, a.site.q	: the number of A sites in the consensus sequence (the aligned region only)
8, t.site.q: the number of T sites in the consensus sequence (the aligned region only)
9, g.site.q: the number of G sites in the consensus sequence (the aligned region only)
10, c.site.q: the number of C sites in the consensus sequence (the aligned region only)
11, ac: the number of A-to-C mutations
12, ag: the number of A-to-G mutations
13, at: the number of A-to-T mutations
14, ca: the number of C-to-A mutations
15, cg: the number of C-to-G mutations
16, ct: the number of C-to-T mutations
17, ga: the number of G-to-A mutations
18, gc: the number of G-to-C mutations
19, gt: the number of G-to-T mutations
20, ta: the number of T-to-A mutations
21, tc: the number of T-to-C mutations
22, tg: the number of T-to-G mutations

"""


import sys,re
argvs = sys.argv

#reverse complement
def inversComplement(input):
    output = ''
    for letter in input:
        letter = letter.lower()

        if letter == 'a':
            output += 't'
        elif letter == 't':
            output += 'a'
        elif letter == 'g':
            output += 'c'
        elif letter == 'c':
            output += 'g'
        else:
            output += letter
    return(output[::-1])

#calculation of atgc composition
def calcNtComposition(l):
  a_num = len([c for c in l if c == 'a'])
  t_num = len([c for c in l if c == 't'])
  g_num = len([c for c in l if c == 'g'])
  c_num = len([c for c in l if c == 'c'])
  compos_l = [a_num,t_num,g_num,c_num]
  return compos_l

#remove gapped/umbigous site
def rmGapSite(s1,s2):
  seq_l = ['a','t','g','c']
  pos_l = [i for i in range(len(s1)) if (s1[i] in seq_l) and (s2[i] in seq_l)]
  s1_nogap_l = [s1[i] for i in range(len(s1)) if i in pos_l]
  s2_nogap_l = [s2[i] for i in range(len(s2)) if i in pos_l]
  return s1_nogap_l,s2_nogap_l

#count Mutations
def countMutations(l1,l2):
  count_d = {'ac': 0, 'gt': 0, 'ag': 0, 'ca': 0, 'cg': 0, 'gc': 0, 'at': 0, 'ga': 0, 'tg': 0, 'ct': 0, 'tc': 0, 'ta': 0}
  mutation_l = ['ac', 'ag', 'at', 'ca', 'cg', 'ct', 'ga', 'gc', 'gt', 'ta', 'tc', 'tg']
  for i in range(len(l1)):
    c1 = l1[i]
    c2 = l2[i]
    if c1 != c2:
      cc = c1 + c2
      count_d[cc] += 1
  count_l = []
  for cc in mutation_l:
    count = count_d[cc]
    count_l.append(count)
  return count_l



#read file
f = open(argvs[1])

align_d = {}
TE_info_l = []
for line in f:
  line = line.rstrip()
  if len (line) > 0:
    initial = line[0]
    initial3 = line[:3]
    if initial.isdigit():
      strand = '+'
      if " C " in line:
        strand = '-'
        line = line.replace(' C ',' ')
      sw,per_div,per_del,per_ins,chrom,start,end,_,TE_name = line.split()[:9]
      Type,Class = TE_name.split('#')
      TE_locus = Class + '|' + Type + '|' + chrom + ':' + start + '-' + end + '|' + strand
      query_name = TE_name
      target_name = chrom
      if len(query_name) > 13:
        query_name = query_name[:13]
      if len(target_name) > 13:
        target_name = target_name[:13]      
      align_d[TE_locus] = {query_name:'',target_name:''}

      TE_info_t = (sw,per_div,per_del,per_ins,chrom,start,end,TE_name,strand,TE_locus,query_name,target_name)
      TE_info_l.append(TE_info_t)

    if (re.match(r'  [^ ]',initial3)) or (re.match(r'C [^ ]',initial3)):
      line = re.sub(r'^C',r'',line)
      line = line.strip()
      line = re.sub('[ ]+',' ',line)
      seq_name,_,seq,_ = line.split(" ")
      align_d[TE_locus][seq_name] = align_d[TE_locus][seq_name] + seq


#if TE insertion is reverse-complement, reverse complement
for TE_locus in align_d:
  strand = TE_locus.split('|')[-1]
  if strand == '-':
    for seq_name in align_d[TE_locus]:
      seq = align_d[TE_locus][seq_name]
      seq = inversComplement(seq)
      align_d[TE_locus][seq_name] = seq



#count, output
mutation_l = ['ac', 'ag', 'at', 'ca', 'cg', 'ct', 'ga', 'gc', 'gt', 'ta', 'tc', 'tg']

print "TE_locus\tchrom\tstart\tend\tsw_score\ttotal_seq_len\ta.site.q\tt.site.q\tg.site.q\tc.site.q\t" + "\t".join(mutation_l) #+ "\tquery_seq\ttarget_seq"

for TE_info_t in TE_info_l:
  sw,per_div,per_del,per_ins,chrom,start,end,TE_name,strand,TE_locus,query_name,target_name = TE_info_t

  query_seq = align_d[TE_locus][query_name].lower()
  target_seq = align_d[TE_locus][target_name].lower()

  query_nogap_l,target_nogap_l = rmGapSite(query_seq,target_seq)
  
  seq_len = len(query_nogap_l)
  
  compos_l = calcNtComposition(query_nogap_l)

  count_l = countMutations(query_nogap_l,target_nogap_l)

  res_l = [TE_locus,chrom,start,end,sw,seq_len] + compos_l + count_l

  print "\t".join([str(c) for c in res_l])

