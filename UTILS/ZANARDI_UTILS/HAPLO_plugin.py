### This program runs utilities for haprep option in 

def hapchk(mapfile,hsnp,cromo):
    ### Read mapfile
    check=open(mapfile).readlines()
    mapp=dict((check[i].strip().split()[1],0) for i in range(len(check)) if check[i].strip().split()[0]==cromo)

    ### Check presence of all SNPs that should be present, otherwise, stop
    missed=[]
    for snp in hsnp:
        if not mapp.has_key(snp):missed.append(snp)
    if missed:return(False,'The following SNPs in BTA19 needed for this analysis are missing:\n            '+'\n            '.join(missed))
    else:return(True,'')

names=[]  ### TEST
def hapconvert(TMPDIR,OUTDIR,output_name,mapfile,pedfile,hsnp,HAPREP):
    haplout=open(OUTDIR+'/HAPREP'+str(output_name)+'.txt','w')
    convert={};alles={};allowed=[];
    haplout.write('ID_animal')
    strandtext={0:'Forward',1:'Top',2:'A/B'}
    for en,a in enumerate(open(mapfile)):
        crom,name,xx,pos=a.strip().split()
        if not hsnp.has_key(name):continue
        names.append(name) #### TEST
        haplout.write(','+name.replace('-','_'))
        convert[(en,0)]={hsnp[name][0]:0,hsnp[name][1]:1}
        convert[(en,1)]={hsnp[name][2]:0,hsnp[name][3]:1}
        convert[(en,2)]={hsnp[name][4]:0,hsnp[name][5]:1}
    haplout.write('\n')
    if HAPREP=='FLK' and en!=1057:
        return(False,'SNPs in '+TMPDIR+'/PLINK_HAPREP.map should be 1058, but got '+str(en+1)+'.\n'+\
             '            Usually this happens if all individauls have missing calls in the "lacking" SNPs.\n'+\
             '            I cannot produce the output file if SNPs are missing, sorry!.')
    strand=99
    for en,a in enumerate(open(pedfile)):
        fid,name,xx,xx,xx,xx,geno=a.strip().split(' ',6)
        geno=geno.strip().split()
        if en==0:   ## Automatic choice of strand
            strands=[0,1,2]
            counts=[0,0,0]
            for x in range(len(geno)/2):
                for ind in strands:
                    try:convert[(x,ind)][geno[x*2]]
                    except: counts[ind]+=1
                    try:convert[(x,ind)][geno[x*2+1]]
                    except: counts[ind]+=1
        for i in range(len(counts)):
            if not counts[i]:
                strand=i
                break
        if strand==99:return(False,"No clear strand was identified. Most probable strand is "+strandtext[counts.index(min(counts))].upper() +"\n"+\
                    "            However, ~"+str(min(counts)/2)+" SNPs are not in the correct strand!\n"+\
                    "            Please check your input file and run again!")
        geno=[str(convert[(x,strand)][geno[x*2]]+convert[(x,strand)][geno[x*2+1]]) for x in range(len(geno)/2)]            
        haplout.write('%s,%s\n' % (name,','.join(geno)))
    haplout.close()
    return(True,strandtext[strand])
