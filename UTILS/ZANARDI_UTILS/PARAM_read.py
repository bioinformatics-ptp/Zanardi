import os

def check_range(inp,min,max,option):
    try:
        val = float(inp)
        min = float(min)
        max = float(max)
    except ValueError:
        return(False,'ERROR in %s: %s is not a number and it should be!' % (option,inp))
    if val<min or val>max:
        return(False,'ERROR in %s: %s is not in the range between %s and %s, but it should be!' % (option,inp,min,max))

##########################
### Checks paths to PLINK
##########################
def check_path(variable,PARAMETERS,DEBUG,exe):
    pgm=variable.strip().split('_')[1]
    matching=[par for par in PARAMETERS if variable in par]
    if not matching:return(False,variable + " required variable not found in parameter file!")
    matching=matching[0].strip().split('=')
    if '~' in matching[1]:return(False,'Zanardi does not accept relative paths (e.g. using ~). Use ONLY absolute paths')
    if len(matching[1]) == 0:
        return(False,pgm+' path not provided in PARAMETER file. If not available, type: python Zanardi.py --download plink')
    else:
        PATH = matching[1].strip()+'/'
        if not os.path.isdir(PATH) or not os.path.exists(PATH+'/'+exe):
            return(False,pgm+' path not provided or executable ( '+exe+' ) not found in: '+PATH+' )')
        else:return(True,PATH)
    if DEBUG:print('PATH SEARCHED FOR PROGMAM '+pgm+' :'+str(PATH))

###############################
### READS INPUT FILE(s) PARAMS
###############################
def inputfi(PARAMETERS,DEBUG,PEDIG):
    pedfile='NOT PROVIDED';mapfile='NOT PROVIDED';Ibull='NOT PROVIDED';Ibull_map='NOT PROVIDED';pedigree='NOT PROVIDED';output_name='';phenotype='NOT PROVIDED'
    variables=['INPUT_PED','INPUT_MAP','INPUT_705','INPUT_705_MAP','INPUT_PEDIG','SPECIES','OUTPUT_NAME','INPUT_PHENO']
    for variable in variables:
        matching=[par for par in PARAMETERS if variable in par][0].strip().split('=')
        if not matching[0]:return(False,variable + " required variable not found in parameter file!")
        if variable == variables[0]:
            if len(matching[1])==0:continue
            pedfile=matching[1].strip().split(',')
            for en,ped in enumerate(pedfile):
                pedfile[en]=ped.strip()
                if not os.path.exists(pedfile[en]):return(False,"File in PLINK fmt (var. INPUT_PED) not found: "+pedfile[en])
        if variable == variables[1]:
            if len(matching[1])==0:continue
            mapfile=matching[1].strip().split(',')
            for en,mapf in enumerate(mapfile):
                mapfile[en]=mapf.strip()
                if not os.path.exists(mapfile[en]):return(False,"Map file in PLINK fmt (var. INPUT_MAP) not found: "+mapfile[en])
        if variable == variables[2]:
            if len(matching[1])==0:continue
            Ibull=matching[1].strip().split(',')
            for en,ib in enumerate(Ibull):
                Ibull[en]=ib.strip()
                if not os.path.exists(Ibull[en]):return(False,"Genotype file in ITB fmt (var. INPUT_705) not found: "+Ibull[en])
        if variable == variables[3]:
            if len(matching[1])==0:continue
            Ibull_map=matching[1].strip().split(',')
            for en,ibmp in enumerate(Ibull_map):
                Ibull_map[en]=ibmp.strip()
                if not os.path.exists(Ibull_map[en]):return(False,"Map file in ITB fmt (var. INPUT_705_MAP) not found: "+Ibull_map[en])
        if variable == variables[4]:
            if len(matching[1])==0 and (PEDIG):return(False,"Pedigree check requested, but no pedigree file provided. Please provide one in paramter file!")
            if len(matching[1])==0:continue
            pedigree=matching[1].strip().split(',')
            if len(pedigree)!=1:return(False,"Only one pedigree file in variable INPUT_PEDIG is accepted!!")
            else:pedigree=pedigree[0]
            if not os.path.exists(pedigree):return(False,"Pedigree file (var. INPUT_PEDIGREE) not found: "+pedigree)
        if variable == variables[7]:
            if len(matching[1])==0:continue
            phenotype=matching[1].strip().split(',')
            if len(phenotype)!=1:return(False,"Only one phenotype file in variable INPUT_PHENO is accepted!!")
            else:phenotype=phenotype[0]
            if not os.path.exists(phenotype):return(False,"Phenotype file (var. INPUT_PHENO) not found: "+phenotype)
        if variable == variables[5]:
            species_dict={'cow':29,'dog':38,'horse':31,'mouse':19,'sheep':26,'human':23,\
                              'pig':18,'goat':29,'all':59,'chicken':38,'wolverine':40}
            species_different=['pig','goat','chicken','human','all']
            if len(matching[1])==0:return(False,"Species variable not found in parameter file!")
            spe=''.join(matching[1].lower().strip().split())
            if not spe in species_dict.keys():
                return(False,"Zanardi does not support this species: %s !!!\nPlease choose one of: %s" %\
                         (spe,' / '.join(species_dict.keys())))
            if not spe in species_different: species=str(' --'+str(spe)+' ')
            if spe in species_different:
                if spe=='human' or spe=='pig' : species=str(' ')
                if spe=='goat' : species=str(' --cow ')
                if spe=='chicken' : species=str(' --dog ')
                if spe=='all' : species=str(' --chr-set 59')
            nchrom=int(species_dict[spe])
        if variable == variables[6]:             #output file name
            matching=[par for par in PARAMETERS if variable in par][0].strip().split('=')
            if len(matching[1])==0:continue
            else:output_name='_'+''.join(matching[1].strip().split())

    if DEBUG:print('INFILE OPTIONS:\n'+'pedfile: '+str(pedfile)+'\n mapfile: '+str(mapfile)+\
                       '\n Ibull: '+str(Ibull)+'\n Ibull_map: '+str(Ibull_map)+'\n pedigree:'+str(pedigree))
    
    return(True,pedfile,mapfile,Ibull,Ibull_map,pedigree,species,nchrom,output_name,phenotype)

### Read parameters for PLINK
def plink_par (PARAMETERS,OUTDIR,output_name,DEBUG):
    variables=['QCMISS_IND','QCMISS_SNP','QCMAF','QCHWE','QC_OTHOPT']
    plinkqc_vals = ['-9','-9','-9','-9','-9','-9']
    for number,variable in enumerate(variables):
        matching=[par for par in PARAMETERS if variable in par][0].strip().split('=')
        if not matching[0]:return(False,variable + " required variable not found in parameter file!")
        if variable in variables[:4]:
            if len(matching[1])==0:continue
            check_range(matching[1],0,1,variable)
            plinkqc_vals[number] = matching[1]
        if variable == variables[4]:             ### Custom (plink) options
            if len(matching[1])!=0:plinkqc_vals[5] = matching[1]
    plinkqc_vals[4] = OUTDIR+'/PLINK_QC'+output_name
    if DEBUG:print('PLINK QC OPTIONS:\n'+str(variables)+'\n'+str(plinkqc_vals))
    return(True,plinkqc_vals)

### Read parameters for PEDIG
def pedig_par(PARAMETERS,DEBUG):
    variables=['PDSKIPCOUPLE','PDMEND_THRES','PDBESTALL']
    for variable in variables:
        matching=[par for par in PARAMETERS if variable in par][0].strip().split('=')
        if not matching[0]:return(False,variable + " required variable not found in parameter file!")
        if variable == variables[0]:
            if len(matching[1])==0:
                skipit = 'NOT PROVIDED'; skipcouples = 0
            else:
                skipcouples = matching[1]
                if not os.path.exists(skipcouples): return(False,"Couples file provided for pedigchk (var. PDSKIPCOUPLE) not found: "+skipcouples)
        if variable == variables[1]:
            checkall= False
            if len(matching[1])==0:return(False,"Mendelian inheritance threshold (var. PDMEND_THRES) not provided! Check PARAMFILE.txt")
            else:
                menderr_thr=matching[1]
                check_range(matching[1],0,1,variable)
        if variable == variables[2]:
            if len(matching[1])!=0:
                if matching[1] == 'Y':
                    ptext = "Checking bestmatch in all individuals in pedigree, as requested"
                    checkall = True
                elif matching[1] == 'N':ptext = "Checking bestmatch in same sex and older than *error* putative parents"
            else: return(False,"Unknown value in parameter file for PDBESTALL variable (expected Y/N, found: "+matching[1]+")!")  
    if DEBUG:print('PEDIG OPTIONS:\n Skipcouples:'+str(skipcouples)+'\n menderr_thr:'+str(menderr_thr))
    return(True,skipcouples,menderr_thr,checkall,ptext)

### Read parameters for MEND
def mend_par(PARAMETERS,DEBUG):
    matching=[par for par in PARAMETERS if 'MENDERR_THRES' in par][0].strip().split('=')
    if not matching[0]:return(False,variable + " required variable not found in parameter file!")
    if len(matching[1])!=0:
        mend_thr=matching[1]
        check_range(matching[1],0,1,'MENDERR_THRES')
    else:return(False,"Mendelian inheritance threshold (var. MENDERR_THRES) not provided! Check PARAMFILE.txt")
    if DEBUG:print('MEND OPTION:\n mend_thres:'+str(mend_thr))
    return(True,mend_thr)

### Read parameters for MDSPLOT
def mds_par(PARAMETERS):
    groupop=False
    matching=[par for par in PARAMETERS if 'MDSGROUPop' in par][0].strip().split('=')
    if not matching[0]:return(False,variable + " required variable not found in parameter file!")
    if len(matching[1])!=0:
        if matching[1] == 'Y':
            mdstext = "Plotting AVERAGE population MDS values"
            groupop = True
        elif matching[1] == 'N':mdstext = "Plotting SINGLE population MDS values"
        else: return(False,"Unknown value in parameter file for MDSGROUPop variable (expected Y/N, found: "+val[1]+")!")  
    return(True,groupop)

### Read beagle v.3 and v.4 params
def beagle_par(PARAMETERS,OUTDIR,output_name,DEBUG):
    variables=['BGMEMORY','BG3_MISSING','BG_OTHOPT']
    beagle_vals= ['2000','0',OUTDIR+'/BEAGLE_OUT','']
    beagle_def=[False,False]

    for variable in variables:
        matching=[par for par in PARAMETERS if variable in par][0].strip().split('=',1)
        if not matching[0]:return(False,variable + " required variable not found in parameter file!")
        if variable in variables[:2]:
            i=variables.index(variable)
            if len(matching[1])==0:continue
            beagle_vals[i]=matching[1]
            beagle_def[i]=True
            check_range(beagle_vals[i],0,9999999999,variable)
        if variable == variables[2]:                 ### Custom (beagle) options
            if len(matching[1])!=0:beagle_vals[3] = matching[1]
    ### Output name file
    beagle_vals[2] = OUTDIR+'/BEAGLE_OUT'+output_name
    if DEBUG:print('BEAGLE OPTIONS:\n'+str(variables)+'\n'+str(beagle_vals))
    return(True,beagle_vals,beagle_def)

### Read FImpute parameters
def fimpute_par(PARAMETERS,OUTDIR,output_name,DEBUG):
    variables=['FMP_NJOB','FMP_OTHOPT']
    fimpute_vals= ['1','']
    fimpute_def=[False,False]

    for variable in variables:
        matching=[par for par in PARAMETERS if variable in par][0].strip().split('=',1)        
        if not matching[0]:
            return(False,variable + " required variable not found in parameter file!")
        if variable in variables[:1]:
            i=variables.index(variable)
            if len(matching[1])==0:continue
            fimpute_vals[i]=matching[1]
            fimpute_def[i]=True
            check_range(fimpute_vals[i],0,9999999999,variable)
        if variable == variables[1]:                 ### Custom (fimpute) options
            if len(matching[1])!=0:fimpute_vals[1] = matching[1]
        
    return(True,fimpute_vals)

### Read ROH parameters
def roh_par(PARAMETERS,TMPDIR,OUTDIR,output_name,DEBUG):
    variables=['ROH_SNP','ROH_MAXMIS','ROH_MAXHET','ROH_MINLEN']
    roh_vals=[15,0,0,1,TMPDIR+'/ROH_reads.txt']

    for variable in variables:
        matching=[par for par in PARAMETERS if variable in par][0].strip().split('=')
        if not matching[0]:return(False,variable + " required variable not found in parameter file!")
        if variable in variables[:4]:
            i=variables.index(variable)
            if len(matching[1])==0:continue
            roh_vals[i]=matching[1]
            check_range(roh_vals[i],0,9999999999,variable)
    roh_vals[4] = OUTDIR+'/ROH_reads'+output_name+'.txt'

    if DEBUG:print('ROH OPTIONS:\n'+str(variables)+'\n'+str(roh_vals))
    return(True,roh_vals)

### Read Admixture vars
def adm_par(PARAMETERS):
    variables=['ADM_KVALUE','ADM_CORE','ADM_CV']
    admixture_vals=[0,1,5]
    for variable in variables:
        matching=[par for par in PARAMETERS if variable in par][0].strip().split('=')
        if not matching[0]:return(False,variable + " required variable not found in parameter file!")
        if variable in variables[1:3]:
            i=variables.index(variable)
            if len(matching[1])==0:continue
            admixture_vals[i]=matching[1]
            check_range(admixture_vals[i],0,9999999999,variable)

        if variable == variables[0]:
            matching=[par for par in PARAMETERS if variable in par][0].strip().split('=')
            if not matching[0]:return(False,variable + " required variable not found in parameter file!")
            if variable in variables[0]:
                i=variables.index(variable)
                if len(matching[1])==0 or int(matching[1]) <= 1 :
                    return(False,variable + " variable not found or minor to 1 in parameter file!")
                admixture_vals[i]=matching[1]
                check_range(admixture_vals[i],0,9999999999,variable)
    return(True,admixture_vals)   
