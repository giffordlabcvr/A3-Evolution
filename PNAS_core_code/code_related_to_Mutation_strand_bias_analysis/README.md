# Analytical codes related to mutation strand bias

## Description
APOBEC3 proteins can suppress the replication of retroviruses including ERVs by inducing G-to-A hypermutations specifically on the positive strand of the viral genome.
To evaluate the degree of APOBEC3 attacks on ERVs or other TEs, these programs evaluate the strand bias of G-to-A mutation rate in each TE subfamily.  
The analytical pipeline is composed of the following three steps:  
  
1) Detecting mutations in respective TE integrants compared to the consensus sequence of the corresponding TE subfamily (done by **detect_mutations_in_each_TE_integrant.py**)  
2) Summarizing mutation information at the TE subfamily level (done by **summarize_mutations_family-level.py**)  
3) Calculating the strand bias of G-to-A mutation rate (done by **calc_mut_rate_family-level.R**)  
  
Detailed information is available in the header of respective programs.
To run these programs on the demo data (**./demo_data/hg19.fa.head_50000.align**), please try the program **batch_demo.sh**.
