# Analytical codes related to the genome screening to extract AID/APOBEC-related sequences
## Description
Two *in-house* programs used in the genome screening (conducted by **Ito et. al., 2019, PNAS**) are available:  
1. select_alignment_site.py  
This program removes the alignment sites in which the site coverage is less than the threshold from multiple sequence alignment.  
Additionally, this program also removes the sequences in which a certain proportion of the alignment sites do not satisfy the criteria of the site coverage.  

2. remove_outlierOTUs.py  
This program removes the OTUs having an extraordinary long external branch in tree from multiple sequence alignment.

Detailed information is available in the header of respective programs.  
To run these programs on the demo data (found in **./demo_data**), please try the program **batch_demo.sh**.  