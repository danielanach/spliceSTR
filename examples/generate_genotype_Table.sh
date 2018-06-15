#!/bin/bash

 
raw=/storage/szfeupe/Runs/650GTEx_estr/Filter_Merged_STRs_All_Samples_New.vcf.gz
Out0='Allele_Gentype_CFTR.table'
Out1='Allele_Gentype_SFTPC.table'
Out2='Allele_Gentype_enos.table'
Out3='Allele_Gentype_ATXN2.table'

# CFTR
python /home/dana/workspace/gtex-estrs/Scripts/Data_preparation/GetNormalizedGenotypes.py --vcf $raw --gtfield 'GB'  --minsamples 200  --region 7:117120017-117308719 --alleles --nonorm  --debug --out $Out0

# SFTPC
python /home/dana/workspace/gtex-estrs/Scripts/Data_preparation/GetNormalizedGenotypes.py --vcf $raw --gtfield 'GB'  --minsamples 200  --region 8:22016184-22021991 --alleles --nonorm --debug  --out $Out1

# enos
python /home/dana/workspace/gtex-estrs/Scripts/Data_preparation/GetNormalizedGenotypes.py --vcf $raw --gtfield 'GB'  --minsamples 200  --region 7:150688144-150711687 --alleles --nonorm --debug   --out $Out2

# ATXN2
python /home/dana/workspace/gtex-estrs/Scripts/Data_preparation/GetNormalizedGenotypes.py --vcf $raw --gtfield 'GB'  --minsamples 200  --region 16:28834356-28848558 --alleles --nonorm --debug  --out $Out3
