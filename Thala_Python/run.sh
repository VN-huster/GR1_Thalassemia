# !/bin/bash

phased=(F01.phased.vcf.gz F02.phased.vcf.gz F07.phased.vcf.gz F08.phased.vcf.gz 
        F09.phased.vcf.gz F10.phased.vcf.gz F12.phased.vcf.gz F13.phased.vcf.gz
        F14.phased.vcf.gz F15.phased.vcf.gz F17.phased.vcf.gz F18.phased.vcf.gz
        F19.phased.vcf.gz F21.phased.vcf.gz F22.phased.vcf.gz F23.phased.vcf.gz
        F24.phased.vcf.gz F25.phased.vcf.gz F26.phased.vcf.gz)

fam=(F01.vcf.gz F02.vcf.gz F07.vcf.gz F08.vcf.gz 
        F09.vcf.gz F10.vcf.gz F12.vcf.gz F13.vcf.gz
        F14.vcf.gz F15.vcf.gz F17.vcf.gz F18.vcf.gz
        F19.vcf.gz F21.vcf.gz F22.vcf.gz F23.vcf.gz
        F24.vcf.gz F25.vcf.gz F26.vcf.gz)

origin_hap=(F01_nipt.pat_origin.hap F02_nipt.pat_origin.hap F07_nipt.pat_origin.hap F08_nipt.pat_origin.hap 
        F09_nipt.pat_origin.hap F10_nipt.pat_origin.hap F12_nipt.pat_origin.hap F13_nipt.pat_origin.hap
        F14_nipt.pat_origin.hap F15_nipt.pat_origin.hap F17_nipt.pat_origin.hap F18_nipt.pat_origin.hap
        F19_nipt.pat_origin.hap F21_nipt.pat_origin.hap F22_nipt.pat_origin.hap F23_nipt.pat_origin.hap
        F24_nipt.pat_origin.hap F25_nipt.pat_origin.hap F26_nipt.pat_origin.hap)

figure=(F01_Sample F02_Sample F07_Sample F08_Sample 
        F09_Sample F10_Sample F12_Sample F13_Sample
        F14_Sample F15_Sample F17_Sample F18_Sample
        F19_Sample F21_Sample F22_Sample F23_Sample
        F24_Sample F25_Sample F26_Sample)

mat_mutation=(5248200 5247153 5247153 5247992 5247992 
            5247992 5248200 5248173 5247992 5247153 
            5247153 5248329 5247992 5248200 5248200 
            5247153 5247992 5247992 5247992)

pat_mutation=(5247992 5248329 5247992 5247153 5247992 
            5247153 5247992 5247992 5248329 5247992 
            5248329 5248200 5247153 5247992 5247992 
            5247153 5248329 5247153 5248329)

len=${#phased[@]}

for((i=0;i<$len;i++))
do
    python3 parent_hap.py "${phased[$i]}" "${fam[$i]}" "${origin_hap[$i]}" "${mat_mutation[$i]}" "${pat_mutation[$i]}"
    python3 famvcf_analysis.py "${phased[$i]}" "${fam[$i]}"
    python3 HMM_pat.py
    python3 HMM_mad.py "${phased[$i]}" "${fam[$i]}" "${origin_hap[$i]}"
    python3 figure.py result/"${figure[$i]}" 
	echo -----------finish "${figure[$i]}"-------------
    
done
