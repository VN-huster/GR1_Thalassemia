with open('vcf.txt','w') as vcf:
    for i in range(10,60):
        vcf.write("bgzip F" + str(i)+".vcf\n")
        vcf.write("tabix F"+str(i)+".vcf.gz\n")