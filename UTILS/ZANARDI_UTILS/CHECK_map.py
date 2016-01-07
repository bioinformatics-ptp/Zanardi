## This program checks the map file, controls consistency and avoids errors related to unknown naming of chromosomes
## G.Marras + E.Nicolazzi

def chk_map(mapfile,nchrom,spe):
    chr_sex={'Y':nchrom+1,'X':nchrom+2,'MT':nchrom+3,'XY':nchrom+4}

    #Read all maps together
    chromosome=dict((map_line.upper().strip().split()[0],0) for map_single in mapfile for map_line in open(map_single))
    war_chr=[True]

    chr_residual = [chr_list for chr_list in chromosome if not chr_list.isdigit()]
    chr_final=list(set(chr_residual)-set(chr_sex))
    if chr_final:return(False,"Chromosome in map file(s) not valid: "+' '.join(chr_final))

#    if [i for i in chr_sex if i in chromosome.keys()]:
    if any(chr_sex.has_key(i) for i in chromosome):
        war_chr= "Provided non-autosomal chromosomes coded as X,Y and/or MT. PLINK will convert them to numbers."
        autosomal=map(int,list(set(chromosome.keys())-set(chr_sex)))
        if nchrom < max(autosomal):
            war_chr= "Number of chromosomes found in map: "+str(max(autosomal))+" | Maximum number of chromosomes for the species chosen: "+str(nchrom)
    else:
#    if not [i for i in chr_sex if i in chromosome]:
        chr_list=sorted([int(k) for k in list(chromosome.keys())])
        if len(chr_list) > nchrom+4 or max(chr_list) > nchrom+4:
            return(False,"Chromosome in map file(s) not valid. Species set: "+spe+" ( "+str(nchrom)+"+4 chr), but found the following: "+','.join(chromosome.keys()))
        if int(chr_list[-1]) in chr_sex.values() :
            war_chr="Non-autosomal chromosomes are coded as numbers. If this is unexpected, please check your map file(s)"
        else:
            war_chr="Non-autosomal chromosomes not found in mapfiles. If this is unexpected, check your map file(s)"
    return(True,war_chr)


def snp_position(mapfile,path):
    error=False;position=[];dict_snp={}
    for line in open(mapfile):
        chrom,SNP,val,pos=line.strip().split()
        if (chrom,pos) in dict_snp:
            if error==False:
                output=open(path+'/Error_SNP_position.txt','w')
                output.write('SNP_ID;CHROM;POSITION\n')
                error=True
            output.write('%s;%s;%s\n' % (SNP,chrom,pos ))
        else:
            dict_snp[chrom,pos]=[0]
    if error:return(False,"Markers with identical CROM+POS found! Imputation software would crash, so I'm stopping.\n"+" "*12+"Please correct (see %s/Error_SNP_position.txt) and rerun " % (path))
    return(True,'')

