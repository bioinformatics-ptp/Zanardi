import subprocess as sub

def allele_freq(PLINK_PATH,pedfile,mapfile,TMPDIR,species,save):
    std = sub.Popen([str(PLINK_PATH+'./plink --noweb '+species+' --ped '+pedfile+' --map '+mapfile+\
          ' --freq --out '+TMPDIR+'/'+save)],shell=True,stdout=sub.PIPE,\
          stderr=sub.STDOUT).stdout.readlines()
    return std


def conversion_PLINK_to_Fimpute(pedfile,save):
    recode={'11':'0','22':'2','12':'1','21':'1','00':'5'}
    output=open(save,'w')
    output.write('ID    Chip                   Call...\n')
    for line in open(pedfile):
        breed,ID,sire,dame,sex,phe,geno=line.strip().split(' ',6)
        geno=geno.split()
        genotype=[recode.get(geno[x]+geno[x+1],'!') for x in range(0,len(geno)-1,2)]
        if genotype.count('!'): return (False,"Error with recode option!! ")
        output.write('%s %s %s\n' % (ID,1,''.join(genotype) ))
    return(True,'ok')

def map_convert_FImipute(mapfile,save):
    output=open(save,'w')
    output.write('SNP_ID Chr Pos Chip1 Chip2\n')    
    conta=1
    for line in open(mapfile):
        chrom,name,val,pos=line.strip().split()
        output.write('%s %s %s %s\n' % (name,chrom,pos,conta))
        conta+=1

def conversion_Fimpute_to_PLINK(imputedfile,TMPDIR,save,pedfile,output_folder,map_name):   
    #geno_ACGT={0:A1A1, 1:A1A2, 2:A2A2, 3:A1A2, 4:A2A1, 5:MISS, 6:MISS, 7:MISS, 8:MISS, 9:MISS}
    geno_ACGT={}
    count=0
    for line in open(TMPDIR+'/freqACGT.frq'):
        if 'MAF' in line:continue
        chrom,SNP,A1,A2,MAF,NCHROBS=line.strip().split()
        if A1=='0': #if snp is monomorphic (A1=0)
            geno_ACGT[SNP]=[A1+' '+A1,A1+' '+A1,A2+' '+A2,A1+' '+A1,A1+' '+A1,A1+' '+A1,A1+' '+A1,A1+' '+A1,A1+' '+A1,A1+' '+A1,SNP]
        else: 
            geno_ACGT[SNP]=[A1+' '+A1,A1+' '+A2,A2+' '+A2,A1+' '+A2,A2+' '+A1,'0 0','0 0','0 0','0 0','0 0',SNP]
        count+=1

    ##keep first 6 column ped file
    ID_info={}
    for line in open(pedfile):
        breed,ID,sire,dame,sex,phe,geno=line.strip().split(' ',6)
        ID_info[ID]=[breed,ID,sire,dame,sex,phe]
    
    #keep snp position
    output_map=open(map_name,'w')
    SNP_info={}
    count=0
    for line in open(output_folder+'/snp_info.txt'):
        if 'SNPID' in line :continue
        SNPID,Chr,BPPos,chip_1=line.strip().split()
        SNP_info[count]=[SNPID]
        count+=1
        output_map.write('%s %s %s %s \n'% (Chr,SNPID,0,BPPos))

    ##conversion FImpute to PLINK format
    output_ped=open(save,'w')
    for line in open(imputedfile):
        lista=[]
        if 'Chip' in line:continue
        sample,chip,geno=line.strip().split()
        for n,a in enumerate(geno):
            lista.append(geno_ACGT[SNP_info[n][0]][int(a)])
        output_ped.write('%s %s \n'% (' '.join(ID_info[sample]),' '.join(lista)))
    return(True)

def pedig_save(pedigree,save):
    output=open(save,'w')
    output.write('ID Sire Dam Sex\n')
    for line in open(pedigree):
        ID,sire,dam,birth,sex=line.strip().split(';')
        output.write('%s %s %s %s\n' % (ID,sire,dam,sex))

def param_FImpute(genotype,snp_info,TMPDIR,output_save,save,OTHOPT):
    output=open(save,'w')
    output.write('title="Zanardi imputation" ;\n') 
    output.write('genotype_file="'+str(genotype)+'" ;\n') 
    output.write('snp_info_file="'+str(snp_info)+'" ;\n')
    output.write('output_folder="'+str(output_save)+'" ;\n')
    output.write('njob='+str(OTHOPT[1][0])+' ;\n') 
    output.write(str(OTHOPT[1][1]).replace(';',';\n')) #OTHER OPTION
