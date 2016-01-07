#program to run pedigree checking over genotypes provided
import sys

## This is the def that actually performs the pedig check
def pedig_check(geno):
    if len(geno[0])!=len(geno[1]):return('ERROR',"Length of genotypes different!. This shouldn't happen!")
    matching=str(sum(geno[0][z]!=geno[1][z] for z in range(len(geno[0])) if geno[0][z]!='.' and geno[1][z]!='.' and geno[0][z]!='-' and geno[1][z]!='-'))
    total=str(sum(geno[0][z]!='.' and geno[1][z]!='.' for z in range(len(geno[0]))))
    return (matching,total) 

def simil_check(geno):
    matching=0;total=0
    for z in range(len(geno[0])):
        if geno[0][z]=='-' or geno[1][z]=='-':continue
        total+=2
        if geno[0][z]!='.' and geno[1][z]!='.':
            if geno[0][z]==geno[1][z]:matching+=2
            continue
        if geno[0][z]=='.' and geno[1][z]=='.':
            matching+=2
            continue
        else:
            matching+=1
    return(matching,total)

def run_mend_check (plink_pedfile, mend_thr, outfile):
    ## Useful dicts 
    rehom={'AA':'A','BB':'B','CC':'C','GG':'G','TT':'T','11':'1','22':'2','00':'-','??':'-','-9-9':'-'}
    mend_thres=float(mend_thr)*100

    ## Output files 
    passped=open(outfile,'w')
    passped.write('%29s; %29s; %12s; %9s; %9s\n' % ('ID_1','ID_2','MendErr/SNPs','MendErr%','Simil%'))

    print "      Doing actual pedigree checks (PASS/FAIL)"
    ## Read genotype file and do the checks - write output files
    cpass=0;checks=0
    iid=[];genotype=[]
    for row in open(plink_pedfile):
        line=row.replace('\t',' ').strip().split()
        idtest=line[0]+' '+line[1]
        genotest=[rehom.get(line[6+x]+line[7+x],'.') for x in range(0,len(line[6:])-1,2)]
        if len(iid)==0:
            genotype.append(genotest)
            iid.append(idtest)
            continue
        for i in range(len(genotype)):
            checks+=1
            if checks%10000==0:print "     %s couples of individuals checked..." % checks
            fast_mend_check=pedig_check([genotest[:2000],genotype[i][:2000]])
            if fast_mend_check[1]=='0':
                print("WARNING: Couple",idtest,"-",iid[i],"cannot be checked. NO homozygote genotypes in common - SKIPPED");continue
            outcome = (float(fast_mend_check[0])/float(fast_mend_check[1]))*100
            if outcome > mend_thres*10.: continue  ### DISCARD CLEARLY WRONG SAMPLES

            ### NOW RUN THE FULL CHECK
            mend_check=pedig_check([genotest,genotype[i]])
            if mend_check[1]=='0':
                print("WARNING: Couple",idtest,"-",iid[i],"cannot be checked. NO homozygote genotypes in common - SKIPPED");continue
            outcome = (float(mend_check[0])/float(mend_check[1]))*100            
            if outcome <= mend_thres:
                cpass+=1
                sim=simil_check([genotest,genotype[i]])
                simil=float(sim[0])/float(sim[1])*100
                passped.write('%29s; %29s; %12s; %9.5f; %9.5f\n' %  (idtest,iid[i], '/'.join(mend_check[:2]),round(outcome,5),round(simil,5)))
        genotype.append(genotest)
        iid.append(line[0]+' '+line[1])

    return (True,str(checks),str(cpass),len(genotype))
