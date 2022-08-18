from cmath import phase
from cyvcf2 import VCF
import sys

def analysis(phased_file, sample_file):
    phased = VCF(phased_file)
    sample = VCF(sample_file)
    dic = {}
    queue = []
    for i in phased:
        queue.append(i.POS)
        dic[i.POS] = [i.genotypes[-1][0], i.genotypes[-1][1], i.genotypes[-2][0], i.genotypes[-2][1] ]
    
    up = 0
    down = 0
    floor_depth = 30
  
    with open('mad.txt','w') as writemad, open('fra.txt','w') as writefra, open('out.txt','w') as writeout, open('ff.txt','w') as writeff:
        writemad.write("CHROM\tPOS\tGenotype\tMat_REF\tMat_ALT\n")
        writefra.write("CHROM\tPOS\tREF\tALT\tFetus_fraction\n")
        writeout.write("CHROM\tPOS\tREF\tALT\tF0\tF1\tM0\tM1\tPlasma_REF\tPlasma_ALT\n")
        for i in sample:
            if (i.POS < queue[0]):
                continue
            elif (i.POS == queue[0]):
                queue.pop(0)

        ## Neu độ sâu bé hơn 30 thì bỏ qua
                if(i.gt_depths[0]<floor_depth or i.gt_depths[1]<floor_depth or i.gt_depths[2]<floor_depth):
                    continue;

        ## Nếu kiểu gen của mẹ dị hợp nhưng tỉ lệ trong khoảng thì bỏ qua
                if(i.genotypes[0]==[0, 1, False]):
                    if( 0.1< min(i.gt_alt_depths[0], i.gt_ref_depths[0])/i.gt_depths[0] < 0.4):  
                        continue;
        ## Nếu kiểu gen của bố dị hợp nhưng tỉ lệ trong khoảng thì bỏ qua
                if(i.genotypes[1]==[0, 1, False]):
                    if( 0.1< min(i.gt_alt_depths[1], i.gt_ref_depths[1])/i.gt_depths[1] < 0.4):
                        continue;

        ## Tính tỉ lệ fetus_fraction
                
            # Mother:Father = 0/0:1/1 
                if(i.genotypes[0]==[0, 0, False] and i.genotypes[1]==[1, 1, False] and 
                (i.genotypes[2]==[0, 0, False] or i.genotypes[2]==[0, 1, False])):
                    fra = 2*i.gt_alt_depths[2]/i.gt_depths[2]
                    writefra.write(i.CHROM+"\t"+str(i.POS)+"\t"+i.REF+"\t"+i.ALT[0]+"\t"+str(fra)+"\n")
                    up += i.gt_alt_depths[2]
                    down += i.gt_depths[2]

            # Mother:Father = 1/1:0/0
                if(i.genotypes[0]==[1, 1, False] and i.genotypes[1]==[0, 0, False] and 
                (i.genotypes[2]==[1, 1, False] or i.genotypes[2]==[0, 1, False]) ):
                    fra = 2*i.gt_ref_depths[2]/i.gt_depths[2]
                    writefra.write(i.CHROM+"\t"+str(i.POS)+"\t"+i.REF+"\t"+i.ALT[0]+"\t"+str(fra)+"\n")
                    up += i.gt_ref_depths[2]
                    down += i.gt_depths[2]
                   
                ## Ghi file Mad
                writemad.write(i.CHROM+"\t"+str(i.POS)+"\t"+str(i.genotypes[0][0])
                        +"/"+str(i.genotypes[0][1])+"\t"+str(i.gt_ref_depths[0])+"\t"+str(i.gt_alt_depths[0])+"\n")

                ## Ghi file out
                writeout.write(i.CHROM+"\t"+str(i.POS)+"\t"+i.REF+"\t"+i.ALT[0]+"\t"+str(dic[i.POS][0])+
                    "\t"+str(dic[i.POS][1])+"\t"+str(dic[i.POS][2])+"\t"+str(dic[i.POS][3])+
                    "\t"+str(i.gt_ref_depths[2])+"\t"+str(i.gt_alt_depths[2])+ "\n")
        ## Ghi file ff
        writeff.write(str(2*up/down))      

phased = sys.argv[1]
fam = sys.argv[2]
analysis(phased, fam)