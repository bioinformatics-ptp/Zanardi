## This is the def that checks PED is sorted correctly (from old to young).
def pedigree_control(pedigree,parents,path,strict):
    diz_ID={'0':'0'};diz_SI_DA={'0':'0'}
    error=False
    Umis=False
    for line in open(pedigree):
        ID,sire,dam,birth,sex=line.strip().split(';')
        if 'UUUUU' in sire:
            sire='0'
            Umis=True
        if 'UUUUU' in dam:
            dam='0'
            Umis=True
        if parents:
            if (sire=='0' and dam != '0') or (dam=='0' and sire != '0'):
                if error==False:
                    output=open(path+'/Error_ID_PEDIGREE.txt','w')
                    output.write('ID_Animal;Sire;Dam\n')
                    error=True
                output.write('%s;%s;%s\n' % (ID,sire,dam ))
        if not diz_ID.has_key(ID):
            diz_ID[ID]=0
            if diz_SI_DA.has_key(ID):return (False,"Pedigree not sorted by date of birth (old to young) - Please sort pedigree (Check animal: "+ID+" )")
            if strict:
                if not diz_ID.has_key(sire):return (False,"Male parent in pedigree file (2nd column) is not present as individual (1st column). Check male parent: "+sire)
                if not diz_ID.has_key(dam):return (False,"Female parent in pedigree file (3rd column) is not present as individual (1st column). Check female parent: "+dam)
                if sex!='M' and sex!='F':return (False,"Unknown sex code in pedigree file (4th column). Individual '"+ID +"' is '"+sex+"' (only 'M' or 'F' are acceptable values)")
            if not diz_SI_DA.has_key(sire):diz_SI_DA[sire]=0
            if not diz_SI_DA.has_key(dam):diz_SI_DA[dam]=0
        else:return(False,"Repeated individual ID in first column - Check sample: "+ID)
    if error:return(False,"BOTH male or female parents (2nd/3rd column) should be present or missing\n"+\
                          " "*12+"Check: %s for full list of individuals with errors." %  (path+'/Error_ID_PEDIGREE.txt'))

    return(True,Umis)

### Renames sire and dam (male and female parent) if ID='UUUUU[...]' (for ITB conversion)
def rename_siredam_pedig(pedigree,path):
    outname=path+'/'+pedigree.strip().split('/')[-1]+'_mod'
    out=open(outname,'w')
    for line in open(pedigree):
        linea=line.strip().split(';')
        if 'UUUUU' in linea[1]:linea[1]='0'
        if 'UUUUU' in linea[2]:linea[2]='0'
        out.write('%s\n' % ';'.join(linea))        
    return(True,outname)

### Checks that all inds genotyped are phenotyped
def id_control_pedig(pedfile,pedig,path):
    error=False
    ID_pedig=dict(((i.strip().split(';')[0]),0) for i in open(pedig)  )
    for line in open(pedfile):
        ID_geno=line.strip().split()[1]
        if not ID_pedig.has_key(ID_geno):
            if error==False:
                output=open(path+'/Error_ID_PEDIGREE.txt','w')
                output.write('ID_not_found_in_pedigree\n')
                error=True
            output.write('%s\n' % (ID_geno))
    return(error,path+'/Error_ID_PEDIGREE.txt')

def id_control_pheno(pedfile,pheno,path):
    error=False
    ID_pheno=dict(((i.strip().split(';')[0]),0) for i in open(pheno)  )

    for line in open(pedfile):
        ID_geno=line.strip().split()[1]
        if not ID_pheno.has_key(ID_geno):
            if error==False:
                output=open(path+'/Error_ID_PHENOTYPES.txt','w')
                output.write('ID_not_found_in_phenotype\n')
                error=True
            else:output.write('%s\n' % (ID_geno))
        return(error,path+'/Error_ID_PHENOTYPES.txt')
