#program to run pedigree checking over genotypes provided
import sys

## This is the def that actually performs the pedig check
def pedig_check(geno):
    if len(geno[0])!=len(geno[1]):return('ERROR',"Length of genotypes different!. This shouldn't happen!")
    matching=str(sum(geno[0][z]!=geno[1][z] for z in range(len(geno[0])) if geno[0][z]!='.' and geno[1][z]!='.' and geno[0][z]!='-' and geno[1][z]!='-'))
    total=str(sum(geno[0][z]!='-' and geno[1][z]!='-' and geno[0][z]!='.' and geno[1][z]!='-' for z in range(len(geno[0]))))
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


def run_pedigree_check (plink_pedfile, pedigree_file, mendel_thres, skipcouples, outfiles, checkall):
    ## Useful dicts 
    rehom={'AA':'A','BB':'B','CC':'C','GG':'G','TT':'T','11':'1','22':'2','00':'-','??':'-','-9-9':'-'}
    sesso={'M':'SIRE','F':'DAM'}
    mendel_thres=float(mendel_thres)*100
    candfound=0

    ## Output files 
    passped=open(outfiles[0],'w')
    failped=open(outfiles[1],'w')
    suggest=open(outfiles[2],'w')
    passped.write('%19s; %19s; %4s; %12s; %9s; %9s\n' % ('ID_son','ID_parent','Type','MendErr/SNPs','MendErr%','Simil%'))
    failped.write('%19s; %19s; %4s; %12s; %9s\n' % ('ID_son','ID_parent','Type','MendErr/SNPs','MendErr%'))
    suggest.write('%12s; %19s; %19s; %4s; %9s; %9s\n' % ('C.Num/T.Cand','ID_son','ID_candidate','Type','MendErr%','Simil%'))

    ## Quick read iids in genotype file (to check all are in pedigree)
    gtyped={}
    for read_ped,row in enumerate(open(plink_pedfile)):
        row=row.replace('\t',' ')
        xx,iid,rest=row.strip().split(' ',2)
        gtyped[iid]=iid

    ## Read the couple checks to avoid (if provided)
    skip_couples={}
    if skipcouples:
        print "      Reading file containing couples to skip..."
        for en,row in enumerate(open(skipcouples)):
            if 'ID_parent' in row:continue
            linea=row.strip().split(';')
            if len(linea)!=5:return(False,'File: '+skipcouples+' is supposed to have 5 (semicolon separated) columns, but found '+str(len(linea)))
            id1,id2,xx,xx,xx=linea
            skip_couples[(id1.strip(),id2.strip())]=0
            skip_couples[(id2.strip(),id1.strip())]=0

    print "      Reading pedigree file and storing matching - genotyped - couples"
    ## Read the pedigree file (just once)
    cfound=0;skipped_comp=0
    couples=[]
    allids={};Yob={};Sob={};found={}
    for read_pedig,row in enumerate(open(pedigree_file)):
        line=row.strip().split(';')
        if len(line)!=5:return(False,'File: '+pedigree_file+' is supposed to have 5 (semicolon separated) columns, but found '+str(len(line)))
        iid=line[0].strip()
        if not gtyped.has_key(iid):continue
        sid=line[1].strip();did=line[2].strip();dob=line[3].strip();sex=line[4].strip()
        if Yob.has_key(iid):return(False,'Individual: '+iid+' found 2 times in pedigree')
        if sex!='M' and sex!='F':return(False,'Error in sex code for individual '+iid+' - M/F allowed, but found '+sex)
        Yob[iid]=int(dob[:4])                        # Year of birth
        Sob[iid]=sex                                 # Sex 
        if gtyped.has_key(sid):                      # If IND+SIRE are genotyped, keep the couple to check
            found[iid]=0
            cfound+=1
            if not skip_couples.has_key((iid,sid)):
                couples.append((iid,sid,'S'))
                found[sid]=0
            else:skipped_comp+=1
        if gtyped.has_key(did):                      # If IND+DAM are genotyped, keep the couple to check
            found[iid]=0
            cfound+=1
            if not skip_couples.has_key((iid,did)):
                couples.append((iid,did,'D'))
                found[did]=0
            else:skipped_comp+=1
            
    print "      Doing actual pedigree checks (PASS/FAIL)"
    ## Read genotype file and do the checks - write output files
    cpass=0;cfail=0;ctrl=0
    search_ped={};pedcheck={}
    orderlist=[]
    for row in open(plink_pedfile):
        line=row.replace('\t',' ').strip().split()
        iid=line[1]
        if not Yob.has_key(iid):return(False,'Individual '+iid+' found in genotype file but not on pedigree file')
        if not found.has_key(iid): continue
        genotype=[rehom.get(line[6+x]+line[7+x],'.') for x in range(0,len(line[6:])-1,2)]      #This zeroes all heterozygous gtypes (reduces calc and storage)
        check_alter=[]
        for i in range(len(couples)):
            if couples[i][0]!=iid and couples[i][1]!=iid:continue
            if pedcheck.has_key(couples[i][0]) and iid==couples[i][1]:
                mend_check=pedig_check([pedcheck[couples[i][0]],genotype])
                if mend_check[1]=='0':
                    print("WARNING: Couple",iid,"-",couples[i][0],"cannot be checked. NO homozygote genotypes in common - SKIPPED");continue
                outcome = (float(mend_check[0])/float(mend_check[1]))*100
                if outcome <= mendel_thres:
                    cpass+=1
                    sim=simil_check([pedcheck[couples[i][0]],genotype])
                    simil=float(sim[0])/float(sim[1])*100
                    if couples[i][2]=='S':passped.write('%19s; %19s; SIRE; %12s; %9.5f; %9.5f\n' %  (couples[i][0],couples[i][1], '/'.join(mend_check[:2]),round(outcome,5),round(simil,5)))
                    else: passped.write('%19s; %19s;  DAM; %12s; %9.5f; %9.5f\n' %  (couples[i][0],couples[i][1], '/'.join(mend_check[:2]),round(outcome,5),round(simil,5)))
                else:
                    cfail+=1
                    if couples[i][2]=='S':
                        failped.write('%19s; %19s; SIRE; %12s; %9.5f\n' %  (couples[i][0],couples[i][1], '/'.join(mend_check[:2]),round(outcome,5)))
                        search_ped[(couples[i][0],'S')]=(pedcheck[couples[i][0]])
                        orderlist.append((couples[i][0],'S'))
                    else:
                        failped.write('%19s; %19s;  DAM; %12s; %9.5f\n' %  (couples[i][0],couples[i][1], '/'.join(mend_check[:2]),round(outcome,5)))
                        search_ped[(couples[i][0],'D')]=(pedcheck[couples[i][0]])
                        orderlist.append((couples[i][0],'D'))
                    continue
            elif pedcheck.has_key(couples[i][1]) and iid==couples[i][0]:
                mend_check=pedig_check([pedcheck[couples[i][1]],genotype])
                if mend_check[1]=='0':
                    print("WARNING: Couple",iid,"-",couples[i][1],"cannot be checked. NO homozygote genotypes in common - SKIPPED");continue
                outcome = (float(mend_check[0])/float(mend_check[1]))*100
                if outcome <= mendel_thres:
                    cpass+=1
                    sim=simil_check([pedcheck[couples[i][1]],genotype])
                    simil=float(sim[0])/float(sim[1])*100
                    if couples[i][2]=='S':passped.write('%19s; %19s; SIRE; %12s; %9.5f; %9.5f\n' %  (couples[i][0],couples[i][1], '/'.join(mend_check[:2]),round(outcome,5),round(simil,5)))
                    else:passped.write('%19s; %19s;  DAM; %12s; %9.5f; %9.5f\n' %  (couples[i][0],couples[i][1], '/'.join(mend_check[:2]),round(outcome,5),round(simil,5)))
                else:
                    cfail+=1
                    if couples[i][2]=='S':
                        failped.write('%19s; %19s; SIRE; %12s; %9.5f\n' %  (couples[i][0],couples[i][1], '/'.join(mend_check[:2]),round(outcome,5)))
                        search_ped[(iid,'S')]=genotype
                        orderlist.append((iid,'S'))
                    else:
                        failped.write('%19s; %19s;  DAM; %12s; %9.5f\n' %  (couples[i][0],couples[i][1], '/'.join(mend_check[:2]),round(outcome,5)))
                        search_ped[(iid,'D')]=genotype
                        orderlist.append((iid,'D'))
                    continue
            else:
                pedcheck[iid]=genotype
                continue

    print "      Searching for best matches on",len(search_ped),"FAILED samples (this step is time consuming depending on # failing samples)\n"
    ## IF there are some failing couples, send them to check procedure
    ## This initializes some useful dicts
    candid_ped={};candid_id={}
    for i in search_ped:
        candid_ped[i]=(100.,0.);candid_id[i]=[]

    ## This does the pedig checking over the available genotypes
    ## Note his is a less efficient -but less memory intensive- way of doing this
    enne=0
    if search_ped:
        for row in open(plink_pedfile):
            line=row.replace('\t',' ').strip().split()
            iid=line[1]
            genotype=[rehom.get(line[6+x]+line[7+x],'.') for x in range(0,len(line[6:])-1,2)]
            for i in search_ped:
                if i[0]==iid:continue                                     #skip same individual check
                if not checkall:
                    if i[1]=='S' and Sob[iid]=='F':continue    #skip dams pop if looking for sires
                    if i[1]=='D' and Sob[iid]=='M':continue    #skip sires pop if looking for dams
                    if Yob[i[0]] <= Yob[iid]:continue                         #skip possible sons
                enne+=1
                ## A first quick check over first 2000 genotypes (looking for 10% error rate)
                fast_check=pedig_check([search_ped[i][:2000],genotype[:2000]])
                if fast_check[1] != '0' and int(fast_check[1]) > 100:  ## If no or low hom in fast comparison, do full comparison (few cases)
                    outcome=(float(fast_check[0])/float(fast_check[1]))*100
                    if outcome > 10.:
                        if outcome < candid_ped[i][0]:candid_ped[i]=(999.,999.);candid_id[i]=[('---',sesso[Sob[iid]])]
                        continue
                ## Full check on passing samples
                itry=pedig_check((search_ped[i],genotype))
                if itry[1]==0:continue
                outcome=float(itry[0])/float(itry[1])*100
                if outcome < candid_ped[i][0]:              #If better match, replace worse match
                    sim=simil_check((search_ped[i],genotype))
                    simil=float(sim[0])/float(sim[1])*100
                    candid_ped[i]=(outcome,simil);candid_id[i]=[(iid,sesso[Sob[iid]])]
                elif outcome == candid_ped[i][0]:           #If same match, add new match 
                    candid_id[i].append((iid,sesso[Sob[iid]]))

        for ind in orderlist:
            tot_cand=str(len(candid_id[ind]))
            for en,cands in enumerate(candid_id[ind]):
                candfound+=1
                num=str(en+1)
                suggest.write('%12s; %19s; %19s; %4s; %9.5f; %9.5f\n' % (num+'/'+tot_cand,ind[0],candid_id[ind][en][0],candid_id[ind][en][1],round(candid_ped[ind][0],5),round(candid_ped[ind][1],5)))
    return (True,str(read_pedig+1),str(read_ped+1),str(len(gtyped)),str(cfound),str(skipped_comp),\
                str(str(cpass+cfail)),str(cpass),str(cfail),candfound,'',enne,len(genotype),enne*len(genotype))
