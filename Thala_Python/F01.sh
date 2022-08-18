#!/bin/bash
python3 parent_hap.py F26.phased.vcf.gz F26.vcf.gz F26_nipt.pat_origin.hap 5247992 5248329
python3 famvcf_analysis.py F26.phased.vcf.gz F26.vcf.gz
python3 HMM_pat.py
python3 HMM_mad.py F26.phased.vcf.gz F26.vcf.gz F26_nipt.pat_origin.hap
python3 figure.py F27_Sample
