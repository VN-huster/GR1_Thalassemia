#Create file haplotype from phased sample
#F01.phased -> parent_haplotype

from cyvcf2 import VCF
import pandas as pd
import sys

def get_haplotype(phased_file, mother_mutation, father_mutation):
    phased = VCF(phased_file)
    h_mat, h_pat = -1, -1
    with open('parent.txt', 'w') as writefile:
        writefile.write("CHROM\tPOS\tREF\tALT\tF0\tF1\tM0\tM1\n")
        for i in phased:
            if i.POS == mother_mutation and i.genotypes[-2][0] == 0:
                h_mat = 1
            if i.POS ==father_mutation and i.genotypes[-1][0] == 0:
                h_pat = 1

            writefile.write(i.CHROM+"\t"+str(i.POS)+"\t"+i.REF+"\t"+i.ALT[0]+"\t"+str(i.genotypes[-1][0])+
          "\t"+str(i.genotypes[-1][1])+"\t"+str(i.genotypes[-2][0])+"\t"+str(i.genotypes[-2][1])+ "\n")
        return (h_mat, h_pat)




def make_haplotype(hap_file, h_mat, h_pat):
    hap = pd.read_csv(hap_file, sep='\t')
    if h_mat == 1 and h_pat == 1:
        hap.rename(columns={'F0':'F1', 'F1':'F0', 'M0':'M1', 'M1':'M0'}, inplace = True)
    elif h_mat ==1 and h_pat == -1:
        hap.rename(columns={'M0':'M1', 'M1':'M0'}, inplace = True)
    elif h_mat ==-1 and h_pat == 1:
        hap.rename(columns={'F0':'F1', 'F1':'F0'}, inplace = True)
    hap.to_csv('Parent.txt', index=False, sep='\t')




phased = sys.argv[1]
arg4 = sys.argv[4]
arg5 = sys.argv[5]
mat, pat = get_haplotype(phased, arg4, arg5)
make_haplotype("parent.txt", mat, pat)