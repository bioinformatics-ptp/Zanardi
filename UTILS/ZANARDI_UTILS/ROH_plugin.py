#Runs of Homozygosity
import subprocess as sub

#####################
### READ MAP FILE ###
#####################
def read_map(map_file,nchrom):
    chromosome=[]
    position=[]
    name=[]
    ##read map file 
    for mapp in open(map_file):
        chro,nam,xx,posiz=mapp.strip().split()
        #if not chro.isdigit() or int(chro)<1 or int(chro)>29: #cambiare i cromosomi!!!
        if not chro.isdigit() or not int(chro) in range(1,nchrom+1): #cambiare i cromosomi!!!
            chromosome.append('skip');position.append('skip')
            continue
        chromosome.append(chro)
        position.append(posiz)
        name.append(nam)
    return(chromosome,position,name)

############################
### RUNS OF HOMOZYGOSITY ###
############################
def ROH(pedfile,list1,list2):

    ##creation of lists
    chromosome,position,name=list1
    mini_SNP1,maxmissing1,maxbuffer1,length_min1,file_out=list2

    ##creation value
    recode={'11':'1','22':'1','12':'2','21':'2','00':'5'}
    breeds={}
    error=False
    mini_SNP=int(mini_SNP1)
    length_min=int(length_min1)
    maxbuffer=int(maxbuffer1)
    maxmissing=int(maxmissing1)

    ##output file
    output=open(file_out,'w')
    output.write('BREED;ANIMAL;CHROMOSOME;COUNT;START;END;LENGTH\n')

    ##def to write to output file
    def write_out(breed,animal,first,last,count,chrom):
        diff=int(last)-int(first)
        if count>mini_SNP and diff/1000000. > length_min:
            output.write('%s;%s;%s;%s;%s;%s;%s\n'%(breed,animal,chrom,count,first,last,diff))
            
    ##identification of ROH
    for en,line in enumerate(open(pedfile)):
        allpp=[]
        breed,ind,sire,dame,sex,phe,geno=line.strip().split(' ',6)  ## NOTE: This only works with SINGLE SPACE - to do (other delims)
        if not breeds.has_key(breed):breeds[breed]=0
        breeds[breed] += 1
        geno=geno.split()
        ##recode genotype 
        genotype=[recode.get(geno[x]+geno[x+1],'!') for x in range(0,len(geno)-1,2)]
        if genotype.count('!'): error=True 
        count1=0;last1=0;first1=0;lastcrom='0'
        for letter in range(len(genotype)):
            if chromosome[letter]=='skip':continue
            if chromosome[letter]!=lastcrom :
                ##first write output
                write_out(breed,ind,first1,last1,count1,lastcrom)
                count1=0
                lastcrom=chromosome[letter]

            if letter >= len(genotype)-1: ##end Chromosome
                if genotype[letter]=='1':
                    write_out(breed,ind,first1,position[letter],count1,lastcrom)
                    continue
                if genotype[letter]=='2' :
                    buff+=1
                    if buff==maxbuffer and missing==maxmissing:
                        write_out(breed,ind,first1,position[letter],count1,lastcrom)
                        continue
                    else:
                        write_out(breed,ind,first1,last1,count1,lastcrom)
                        continue
                if  genotype[letter]=='5':
                    missing+=1
                    if buff==maxbuffer and missing==maxmissing:
                        write_out(breed,ind,first1,position[letter],count1,lastcrom)
                        continue
                    else:
                        write_out(breed,ind,first1,last1,count1,lastcrom)
                        continue

            if count1==0:
                buff=0;first1=0; last1=0; missing=0;
                if genotype[letter]=='1':
                    count1+=1
                    first1=position[letter]
                    continue
                else:continue
            if genotype[letter]=='1':
                count1+=1
                last1=position[letter]
                continue

            elif genotype[letter]=='2': ##Control for heterozygous genotypes
                buff+=1
                if buff<=maxbuffer and missing<=maxmissing:
                    count1+=1;continue
                elif buff > maxbuffer:
                    last1=position[letter-1]
                    ##second write output
                    write_out(breed,ind,first1,last1,count1,lastcrom)
                    count1=0
                    continue
            elif genotype[letter]=='5': ##Control for missings genotypes
                missing+=1
                if buff<=maxbuffer and missing<=maxmissing:
                    count1+=1;continue
                elif missing > maxmissing:
                    last1=position[letter-1]
                    ##third write output
                    write_out(breed,ind,first1,last1,count1,lastcrom)
                    count1=0

    return (breeds,error)

##############################
### INBREEDING COEFFICIENT ###
##############################
def inbreeding_ROH(runs_file,map_file,save,nchrom):
    error=False
    length={}
    n_crom='1'
    ##read map file
    for line in open(map_file):
        chromo,nam,xx,posiz=line.strip().split()
        if not chromo.isdigit():continue
        if not int(chromo) in range(1,nchrom+1,1):continue
        if n_crom==int(chromo):length[chromo]=int(posiz)
        n_crom=int(chromo)

    genome=sum(length.values())
    animal_sumroh={};animal_chrroh={}
    for line in open(runs_file):
        if 'BREED' in line:continue
        breed,id_animal,chromosome,count,start,end,difference=line.strip().split(';')
        if not animal_sumroh.has_key((breed,id_animal)):animal_sumroh[(breed,id_animal)]=0
        if not animal_chrroh.has_key((breed,id_animal,chromosome)):animal_chrroh[(breed,id_animal,chromosome)]=0
        animal_sumroh[(breed,id_animal)]+=float(difference)
        animal_chrroh[(breed,id_animal,chromosome)]+=float(difference)
        
    output=open(save,'w')
    allchrom=sorted([int(x) for x in length.keys()])
    output.write('BREED;ANIMAL_ID;FROH_ALL;FROH_CHR%s \n' % ';FROH_CHR'.join(str(x) for x in allchrom))
    for runs in animal_sumroh: 
        FbyChrom=[str(round(animal_chrroh.get((runs[0],runs[1],str(x)),0)/float(length[str(x)]),4)) for x in allchrom]       
        output.write('%s;%s;%s\n' % (';'.join(runs),animal_sumroh[runs]/genome, ';'.join(FbyChrom)))
    
    
############################    
### COUNT SNP INSIDE ROH ###
############################
def snp_inside_ROH(dati_roh,map_file,save,n_breed,nchrom):

    save=open(save,'w')
    save.write('SNP_NAME;CHR;POSITION;COUNT;BREED;PERCENTAGE\n')

    animal={};marker={};final_SNP={};snp_name={}
    for val in open(map_file,'r'):
        chrom,name,x,pos=val.strip().split()
        if int(chrom) > nchrom or int(chrom) < 1:continue 
        if not marker.has_key(chrom):
            marker[chrom]=[]
            final_SNP[chrom]=[]
            snp_name[chrom]=[]
        marker[chrom].append(int(pos))
        final_SNP[chrom].append(0)
        snp_name[chrom].append(name)

    for bre in n_breed.keys():
        n_animal=int(n_breed[bre])
        for read in open(dati_roh,'r'):
            if 'CHROMOSOME' in read:continue
            breed,animal,chrom,conta,inizio,fine,differenza=read.strip().split(';')

####?!?!?!?!? GABRIELE CONTROLLA ##########
            if chrom=='30' or breed !=bre:continue ### THIS WORKS ONLY ON COW
            #if int(chrom)<=nchrom or breed !=bre:continue ### THIS WORKS ALWAYS
###########################################

            start=int(inizio);end=int(fine)
            for val,pos in enumerate(marker[chrom]):
                if pos < start:continue
                if pos > end:break
                final_SNP[chrom][val] +=1    

        for t in final_SNP:
            for n,y in enumerate(final_SNP[t]):
                save.write('%s;%s;%s;%s;%s;%s\n' % (snp_name[t][n],t,marker[t][n],final_SNP[t][n],bre,int(final_SNP[t][n])*100/int(n_animal)))
                final_SNP[t][n]=0


#############################
#### PLOT SNP INSIDE ROH ####
#############################
#def runMDS_R(TMPDIR,OUTDIR,groupop)
def plot_ROH(TMPDIR,OUTDIR,name_file):
    OUTroh=TMPDIR+"/SNP_inside_roh.txt"
    roh=open(TMPDIR+'/roh_plot.R','w')
    roh.write("setwd('"+TMPDIR+"')\n")
    roh.write("library('ggplot2')\n")
    roh.write("all_breed<-read.table('"+OUTroh+"',header=T,sep=';')\n")
    roh.write("pdf('"+OUTDIR+"/ROH_plot"+name_file+".pdf',height=12, width=20)\n")    
    roh.write("head(all_breed) \n")
    roh.write("seq(as.factor(table(all_breed$CHR))) \n")
    roh.write("for (a in levels(as.factor(all_breed$CHR))){  \n cromo<-subset(all_breed,CHR==a) \n")
    roh.write("grafico=ggplot(data=cromo,aes(x=POSITION,y=PERCENTAGE,colour=BREED))  + geom_line() +  ggtitle(paste('chr',a,sep=' ')) + scale_y_continuous(limits = c(-0, 100)) \n")
    roh.write("print(grafico) } \n dev.off() \n")
    roh.write("print('ENDOK')\n")
    roh.close()
    
    ##LAUNCH PGM
    sub.call(["R CMD BATCH %s/roh_plot.R %s/roh_plot.Rout " % (TMPDIR,TMPDIR)],shell=True)

    return (str("/ROH_plot"+name_file+".pdf"))
