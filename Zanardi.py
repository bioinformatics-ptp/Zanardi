"""
Program that easily launches genetic programs in a modular fashion.
Original version: Ezequiel L. Nicolazzi & Gabriele Marras (Parco Tecnologico Padano, Italy), 2015.

THIS PROGRAM RUNS WITH PYTHON >2.5. ANY VERSION < 2.5 WILL NOT WORK!!!

NOTE FOR PYTHONISTS:
      1) Most of internal checks could be implemented using the native "assert" command.
         We chose not to use this HIGHLY useful command to add some more info (in a neat way) for non-python users.
      2) We're aware there are several libraries able to perform many of these tasks in a more efficient fashion (e.g. numpy).
         The idea is to avoid all non-native python libraries, in order to avoid users installing stuff they might not use..
      3) Apart from that, more stylish and efficient programming comments/advices/tips are always welcome!!

History:
 - 20150408: GM+ELN: First alpha-release of the program 
   20150508: GM+ELN: Adding --pedigchk, --roh and --froh options
   20150612: GM+ELN: Adding --admixture and chromosome check
   20150625: GM+ELN: Adding --haprep option
   20150829: ELN+GM: Change the structure. Module all activities and minor bug fix.
   20150925: ELN: Improve speed in option --pedigchk
   20151006: GM+ELN: Bug fix (beagle3 and beagle4)
   20151008: GM : Adding --fimpute option
   20151209: ELN: Adding --mend option 
   20170623: GM: Bug fix (--roh), Update PLINK link

For bug report/comments, please write to: ezequiel.nicolazzi@ptp.it or gabriele.marras@ptp.it
"""

### Import python modules
import sys,os,time,gzip
from optparse import OptionParser
import subprocess as sub

### Create new path for own modules
pydir_util = os.path.abspath(os.path.join(os.path.dirname(__file__),"./UTILS/ZANARDI_UTILS"))
if not pydir_util in sys.path:sys.path.insert(1,pydir_util)

### Import own modules
import DOWNLOAD_plugin, PARAM_read, CHECK_map,\
       PLINK_plugin, INTERBULL_convert, HAPLO_plugin,\
       PEDIG_plugin, ROH_plugin, ADMIXTURE_plugin,\
       PHENO_PEDIG_plugin, FIMPUTE_plugin, MEND_plugin

################################################
### Useful defs
################################################
# Prints (if required) and writes the log 
def logit(msg):
    global VER
    if VER:print(msg)
    log.write(msg+'\n')

# This stops the pipeline when something BAD happens
def bomb(msg):
    logit('\n\n[BAD NEWS]: '+msg+'\n')
    logit("*"*60+'\nProgram stopped because of an error. Please read .log file.\n'+"*"*60)
    sys.exit()

def warn(msg):
    logit('     =======')
    logit(' ==> WARNING <==: %s' % msg)
    logit('     =======      ')    

################################################
### Options parser 
################################################
analysis=False
HERE = os.getcwd()
UTILDIR = HERE+'/UTILS/'
usage = "usage: %prog [options]"
parser = OptionParser(usage=usage)

parser.add_option("--param",dest = "PARAM",default = HERE+'/PARAMFILE.txt',
                  help = "Changes name/location of default param file")
parser.add_option("--download",dest = "DOWNLOAD", default = '',
                  help = "Download programs required for the analyses - See README file")
parser.add_option("--plinkqc", action = "store_true", dest = "PLINKQC", default = False,
                  help = "Runs PLINK to QC your data")
parser.add_option("--mds", action = "store_true", dest = "MDSPLOT", default = False,
                  help = "Runs MDS analysis using PLINK and R (lib. ggplot2)")
parser.add_option("--pedigchk", action = "store_true", dest = "PEDIG", default = False,
                  help = "Runs pedigree checks (using M.Errors)")
parser.add_option("--mendchk", action = "store_true", dest = "MEND", default = False,
                  help = "Runs a mendelian inheritance check over all gtyped samples")
parser.add_option("--beagle3", action = "store_true", dest = "BEAGLE3", default = False,
                  help = "Runs BEAGLE v.3 imputation software")
parser.add_option("--beagle4", action = "store_true", dest = "BEAGLE4", default = False,
                  help = "Runs BEAGLE v.4 imputation software")
parser.add_option("--fimpute", action = "store_true", dest = "FIMPUTE", default = False,
                  help = "Runs FImpute imputation software")
parser.add_option("--haprep", dest = "HAPREP", default = False,
                  help = "Prepares files for haplotype check of carriers - provide breed (FLK or BSW)")
parser.add_option("--froh", action = "store_true", dest = "FROH", default = False,
                  help = "Coefficient Inbreeding from ROH")
parser.add_option("--roh", action = "store_true", dest = "ROH", default = False,
                  help = "Runs of Homozygosity")
parser.add_option("--admixture", action = "store_true", dest = "ADMIXTURE", default = False,
                  help = "Runs ADMIXTURE software ")
parser.add_option("--save",action = "store_true", dest = "SAVEIT", default = False,
                  help = "Avoids deleting TEMP folder")
parser.add_option("--outdir",dest = "OUTDIR",default = HERE+"/OUTPUT",
                  help = "Change default name/location OUTPUT directory")
parser.add_option("--tmpdir",dest = "TMPDIR",default = HERE+"/TEMP",
                  help = "Change default name/location TEMP directory")
parser.add_option("--debug", action = "store_true", dest = "DEBUG", default = False,
                  help = "Use this for debugging" )
parser.add_option("-q","--quiet",action = "store_false", dest = "VER", default = True,
                  help = "Avoid showing runtime messages to stdout")

(opt, args) = parser.parse_args()
VER = opt.VER
log = open(sys.argv[0].strip().split('.')[0]+'.log','w')

_ALLOPTS_=(opt.PLINKQC,opt.MDSPLOT,opt.PEDIG,opt.MEND,opt.BEAGLE3,opt.BEAGLE4,opt.FROH,opt.ROH,opt.ADMIXTURE,opt.HAPREP,opt.FIMPUTE)
_IMPUTATION_=(opt.BEAGLE3,opt.BEAGLE4,opt.FIMPUTE)

##########################################
### Main (and raw) system check
##########################################
# Extract environment info (system and memory) to try avoiding some sys-dpd errors. ############
SyS = sub.Popen(["uname -s"], shell=True,stdout=sub.PIPE, stderr=sub.STDOUT).stdout.readline().strip()
if SyS =='Darwin':
    Sys = 'Mac';
    logit('Recognised a MAC OSx system')
    warn('MAC OSx support is still BETA!')
else: 
    Sys = 'Linux';logit('Recognised a Linux/Linux-like system')

###########################################################################
### If --download option is required, do this first and then exit program
##########################################################################
if opt.DOWNLOAD:
    usropt=opt.DOWNLOAD.lower().strip().split(',')
    logit("Please wait... downloading "+" and ".join(usropt))
    down_stuff=DOWNLOAD_plugin.runit(usropt,UTILDIR,SyS,HERE,opt.DEBUG)
    if down_stuff[0]:
        logit(down_stuff[1]);sys.exit()
    else:bomb(down_stuff[1])

###################################################################################
### Create temporary and output folders if not present, check some essential stuff
###################################################################################
sub.call(["mkdir -p "+opt.TMPDIR],shell = True)
sub.call(["mkdir -p %s" % opt.OUTDIR],shell = True)
if not os.path.isdir(opt.TMPDIR):bomb("Oops! Couldn't create TEMP directory")
if not os.path.isdir(opt.OUTDIR):bomb("Oops! Couldn't create output directory ( %s )" % opt.OUTDIR)
if not os.path.exists(opt.PARAM):bomb("Parameter file not found!!! Looking for: %s" % opt.PARAM)

if not (_ALLOPTS_):bomb("An analysis is required! Please type 'python Zanardi.py -h' for a list of options")
if sum(_IMPUTATION_)>1:bomb("Imputation procedures (Beagle3/4 + FImpute) cannot be run simultaneously!")

### Check haprep options. And if no imputation routine is chosen, activate BEAGLE4 by default.
if opt.HAPREP:
    if sum(_IMPUTATION_)!= 0:bomb("Option --haprep does not accept imputation software option")
    if not opt.HAPREP in ('BSW','FLK'): bomb('Option --haprep accepts only BSW (Brown Swiss) and FLK (Fleckvieh) parameters. NOTE: Use only these breeds!')

########################################################
# READING AND IMPORTING PARAMFILE options (and check it)
########################################################
PARAMFILE_OPEN=open(opt.PARAM,'r').readlines()
PARAMETERS=[ind.strip().replace("'","").replace('"','').replace(';',',') for ind in PARAMFILE_OPEN if not '#' in ind[0]]
if opt.DEBUG:logit('PARAMETER FILE:'+opt.PARAM+':\n'+str(PARAMETERS))

########################################################
# A) PARAM READ: CHECKING PATHS TO 3rd PARTY SOFTWARE
########################################################

### PLINK (required for all analyses)
cpath=PARAM_read.check_path('PGM_PLINK',PARAMETERS,opt.DEBUG,'plink')
if cpath[0]:PLINK_PATH=cpath[1]
else:bomb(cpath[1])

### Third party software check
if opt.BEAGLE3:
    cpath=PARAM_read.check_path('PGM_FCGENE',PARAMETERS,opt.DEBUG,'fcgene')
    if cpath[0]:FCGENE_PATH=cpath[1]
    else:bomb(cpath[1])
    cpath=PARAM_read.check_path('PGM_BEAGLE3',PARAMETERS,opt.DEBUG,'beagle3.jar')
    if cpath[0]:BEAGLE_PATH=cpath[1]
    else:bomb(cpath[1])

if opt.BEAGLE4:
    cpath=PARAM_read.check_path('PGM_BEAGLE4',PARAMETERS,opt.DEBUG,'beagle4.jar')
    if cpath[0]:BEAGLE_PATH=cpath[1]
    else:bomb(cpath[1])

if opt.ADMIXTURE:
    cpath=PARAM_read.check_path('PGM_ADMIXTURE',PARAMETERS,opt.DEBUG,'admixture')
    if cpath[0]:ADMIXTURE_PATH=cpath[1]
    else:bomb(cpath[1])

if opt.FIMPUTE:
    cpath=PARAM_read.check_path('PGM_FIMPUTE',PARAMETERS,opt.DEBUG,'FImpute')
    if cpath[0]:FIMPUTE_PATH=cpath[1]
    else:bomb(cpath[1])

########################################################
# B) PARAM READ: OPTION PARAMETERS READ
########################################################

### PARAM READ: INPUT FILES
inputf=PARAM_read.inputfi(PARAMETERS,opt.DEBUG,opt.PEDIG) 
if inputf[0]:pedfile,mapfile,Ibull,Ibull_map,pedigree,spe,nchrom,output_name,phenotype=inputf[1:]
else:bomb(inputf[1])

### PARAM READ: OPTS FOR PLINK QC
if opt.PLINKQC: 
    ppar=PARAM_read.plink_par(PARAMETERS,opt.OUTDIR,output_name,opt.DEBUG)
    if ppar[0]:plinkqc_vals=ppar[1]
    else:bomb(ppar[1])

### PARAM READ: OPTS FOR PEDIGCHK
if opt.PEDIG: 
    pedpar=PARAM_read.pedig_par(PARAMETERS,opt.DEBUG)
    if pedpar[0]:skipcouples,menderr_thr,checkall,ptext=pedpar[1:]
    else:bomb(pedpar[1])

### PARAM READ: OPTS FOR MENDCHK
if opt.MEND: 
    mendpar=PARAM_read.mend_par(PARAMETERS,opt.DEBUG)
    if mendpar[0]:mend_thr=mendpar[1]
    else:bomb(mendpar[1])

### PARAM READ: OPTS FOR MDS PLOT
if opt.MDSPLOT: 
    mdspar=PARAM_read.mds_par(PARAMETERS)
    if mdspar[0]:groupop=mdspar[1]
    else:bomb(mdspar[1])

### PARAM READ: OPTS FOR BEAGLE 3 and 4
if opt.BEAGLE3 or opt.BEAGLE4: 
    beapar=PARAM_read.beagle_par(PARAMETERS,opt.OUTDIR,output_name,opt.DEBUG)
    if beapar[0]:beagle_vals,beagle_def=beapar[1:]
    else:bomb(beapar[1])

### PARAM READ: OPTS FOR FIMPUTE
if opt.FIMPUTE: 
    PARAMETERS_fimpute=[ind.strip() for ind in PARAMFILE_OPEN if not '#' in ind[0]]
    fimpar=PARAM_read.fimpute_par(PARAMETERS_fimpute,opt.OUTDIR,output_name,opt.DEBUG)
    if fimpar[0]:fimpute_vals,fimpute_def=fimpar[1]
    else:bomb(fimpar[1])

### PARAM READ: OPTS FOR ROH & FROH
if opt.FROH or opt.ROH: 
    rohpar=PARAM_read.roh_par(PARAMETERS,opt.TMPDIR,opt.OUTDIR,output_name,opt.DEBUG)
    if rohpar[0]:roh_vals=rohpar[1]
    else:bomb(rohpar[1])

### PARAM READ: OPTS FOR ADMIXTURE
if opt.ADMIXTURE: 
    admpar=PARAM_read.adm_par(PARAMETERS)
    if admpar[0]:admixture_vals=admpar[1]
    else:bomb(admpar[1])

#######################################################
# Print and check options before running the program
#######################################################
logit('\n'+'*'*81+'\n'+'*'*81)
logit('*'*32+'  Zanardi Suite  '+'*'*32)
logit('*'*28+'  FP7 Gene2farm project  '+'*'*28)
logit('*'*81+'\n'+'*'*81)
logit(' ==> PROGRAM STARTS: '+time.strftime("%B %d, %Y - %l:%M%p %Z")+'\n'+'-'*81)
logit('\n ==>  %-30s  %-s' % ('Parameter file used            :',opt.PARAM))
if opt.PLINKQC or opt.MDSPLOT or opt.BEAGLE3:
    logit(' ==>  %-30s  %-s' % ('Path to PLINK pgm provided     :', PLINK_PATH))
if opt.BEAGLE3:
    logit(' ==>  %-30s  %-s' % ('Path to fcGENE pgm provided    :', FCGENE_PATH))
    logit(' ==>  %-30s  %-s' % ('Path to BEAGLE v.3 pgm provided:', BEAGLE_PATH))
if opt.BEAGLE4:logit(' ==>  %-30s   %-s' % ( 'Path to BEAGLE v.4 pgm provided:', BEAGLE_PATH))
logit(' ==>  %-30s   %-s' % ('Species selected for PLINK     :',spe.upper().replace('--','').strip()))
logit(' ==>  %-30s' % 'Input files provided')
logit('      - %-30s %-s' % ( 'PED file                     :',str(pedfile).replace('[','').replace(']','')))
logit('      - %-30s %-s' % ( 'MAP file                     :',str(mapfile).replace('[','').replace(']','')))
logit('      - %-30s %-s' % ( 'ITB geno file                :',str(Ibull).replace('[','').replace(']','')))
logit('      - %-30s %-s' % ( 'ITB map file                 :',str(Ibull_map).replace('[','').replace(']','')))
logit("      - %-30s '%-s'" % ( 'PEDIGREE file                :',str(pedigree)))
logit("      - %-30s '%-s'" % ( 'PHENOTYPE file               :',str(phenotype)))
if output_name:logit('      - %-30s %-s\n' % ( 'OUTPUT prefix                :',output_name.replace('_','')))
else:logit('      - %-30s %-s\n' % ( 'OUTPUT prefix                :','NOT PROVIDED'))

if opt.MDSPLOT or opt.ROH or opt.ADMIXTURE: 
    logit('      NOTE: Options --mds, --roh and --admixture require ggplot2 R library installed. If not installed please run,')
    logit('            from within R framework, the following command (having root privileges):')
    logit("             -> install.packages('ggplot2')")
    logit("")

if spe.upper().replace('--','').strip() != 'CHR-SET 59':
    if Ibull != "NOT PROVIDED"  and spe.upper().replace('--','').strip() != 'COW':
        bomb('if you provide 705 file, change SPECIES in COW or ALL in parameter file!!')

####################################################
### DOING STUFF: Interbull conversion (if req'd)
####################################################
if pedfile == "NOT PROVIDED":pedfile=[];mapfile=[]
if Ibull != "NOT PROVIDED":  
    logit(" ==>  Converting Interbull 705 fmt into PLINK (PED & MAP) fmt \n")
    tot_ani=0

    ## Do the job
    outcome = INTERBULL_convert.i705_convert_to_plink(Ibull,Ibull_map,opt.TMPDIR+'/convert_705')
    if opt.DEBUG:logit('OUTCOME_INTERBULL CONVERSION:\n'+str(outcome)+'\n')

    ## Report possible error (map with spaces will make plink crash) 
    warn_maps=outcome[3]
    for chips in warn_maps:
        if warn_maps[chips]:bomb("SNP IDs with blank spaces found in SNP map (not allowed):"+str(chips))

    ## Report possible error (Individuals discarded for issues with the number of SNPs)
    if not outcome[0]:
        name_errorfile='/ITB_to_PLINK_unconverted.txt'
        warn("The following individuals were excluded during the conversion process (not sorted) - see "+opt.OUTDIR+name_errorfile+":")
        e_rep=outcome[4]
        for indiv in e_rep:
            if isinstance(e_rep[indiv][2],bool):
                if e_rep[indiv][2]:logit(' '*11+'ID: '+indiv+' - No 705_MAP corresponding to # SNPs in genotype count ( '+e_rep[indiv][0]+')')
                else:logit(' '*11+'ID: '+indiv+' - Wrong number of SNPs in 705 file. Expected: '+e_rep[indiv][0]+' - Found: '+e_rep[indiv][1])
            else:
                bomb(' '*11+'ID:'+indiv+' - Wrong sex in 705 file. Expected: M/F/U - Found: '+e_rep[indiv][1])
        save=open(opt.OUTDIR+name_errorfile,'w')
        save.write( 'ANIMAL ID \n %s' %  ('\n'.join(e_rep.keys())))
        save.close()
        
    ## Report converted individuals
    for n_convert in outcome[2]:
        tot_ani+=outcome[2][n_convert]
        logit(" "*26+"Map size: %s SNPs -- Animals converted: %s " % (n_convert,outcome[2][n_convert]) )
    logit("\n"+" "*14+"[GOOD NEWS]: 705 file converted!! \n"+" "*26 +" Total animals converted: %s\n" % tot_ani)

    for chips in outcome[1]:
        if os.path.exists(opt.TMPDIR+'/convert_705_'+chips+'.ped'):
            pedfile.append(opt.TMPDIR+'/convert_705_'+chips+'.ped')
            mapfile.append(opt.TMPDIR+'/convert_705_'+chips+'.map')
        else:
            bomb('For some reason, file: '+opt.TMPDIR+save_705+'_'+chips+'.ped was not created. Please check!')
    if opt.DEBUG:logit('PED_MAP_AFTER_INTERBULL:\n Pedfile:\n'+str(pedfile)+'\n Mapfile:\n'+str(mapfile))


### Check mapfile and species
mapchk=CHECK_map.chk_map(mapfile,nchrom,spe) 
if not mapchk[0]:bomb(mapchk[1])
if mapchk[1]:warn(mapchk[1])

#####################################################################################################
### DOING STUFF: PLINK merge of all plinks available (if req'd) and standardize separator (if req'd)
#####################################################################################################
mergeit=False
if len(pedfile)==1: 
    pedfile=pedfile[0];mapfile=mapfile[0]
    if len(open(pedfile).readline().strip().split('\t')) > 6:  #If separator is not blank space
        warn("Separator is not blank space. Now standardizing separator...")
        outcome = PLINK_plugin.runSTD_PLINK(PLINK_PATH,pedfile,mapfile,opt.TMPDIR,spe,opt.OUTDIR)
        if any('Error:' in string for string in outcome):
            bomb("PLINK STD failed! Please check %s/PLINK_STND.log to find what went wrong!" % (opt.OUTDIR))
        else:
            logit("\n"+" "*14+"[GOOD NEWS]: PLINK STANDARDIZATION successful!")
            logit(" "*27+              "Output PLINK (standardized)     : "+opt.OUTDIR+"/PLINK_STND[.map/.ped] \n")
            pedfile=opt.OUTDIR+'/PLINK_STND.ped';mapfile=opt.OUTDIR+'/PLINK_STND.map'
else:
    mergeit=True
    if len(pedfile)!=len(mapfile): bomb('The number of PED and MAP files provided must be the same!')
    warn('Multiple genotype files provided, thus default --merge PLINK option will be used.\n'+\
          ' '*18+'Only CONSENSUS calls will be retained (e.g. setting to missing all non consensus)')
if opt.DEBUG:logit('PED_MAP_BEFORE_MERGE:\n Pedfile:\n'+str(pedfile)+'\n Mapfile:\n'+str(mapfile))

### If need to merge, do the job
if mergeit:
    outcome = PLINK_plugin.runMERGE_PLINK(PLINK_PATH,pedfile,mapfile,opt.TMPDIR,spe,opt.OUTDIR)
    ## Checks and error report
    if opt.DEBUG:logit('OUTCOME_PLINK_MERGE:\n'+str(outcome)+'\n')
    if any('Error:' in string for string in outcome):
        bomb("PLINK MERGE failed! Please check %s/PLINK_MERGED.log to find what went wrong!" % (opt.OUTDIR))
    else:
        logit("\n"+" "*14+"[GOOD NEWS]: PLINK merge successful!")
        logit(" "*31+              "- Output PLINK (merged)     : "+opt.OUTDIR+"/PLINK_MERGED[.map/.ped] \n")
        pedfile=opt.OUTDIR+'/PLINK_MERGED.ped';mapfile=opt.OUTDIR+'/PLINK_MERGED.map'   
    if any('Warning' in string for string in outcome):
        warn("Presence of markers with identical position/chrom found.. if this is unexpected STOP the program and check!")

if opt.DEBUG:logit('PED_MAP_AFTER_MERGE:\n Pedfile:\n'+str(pedfile)+'\n Mapfile:\n'+str(mapfile))

###################################################################################
### DOING STUFF: Check error in pedrigree/phenotype file for PEDIG option 
###################################################################################
## PEDIGREE check
fimpute_pedigree=False
if opt.FIMPUTE and pedigree!='NOT PROVIDED':fimpute_pedigree=True

if opt.PEDIG or fimpute_pedigree:
    logit('\n      -------------------')
    logit(' ***  PEDIGREE FILE CHECK ***')
    logit('      -------------------')
    #control ID in geno-pedigree file
    outcome = PHENO_PEDIG_plugin.id_control_pedig(pedfile,pedigree,opt.OUTDIR)
    if outcome[0]:
        bomb("Some genotyped individuals found missing in pedigree file.\n"+" "*\
                 12+"Check full list in %s" % ( outcome[1]) )
    #control sort and consistency
    if opt.PEDIG: outcome = PHENO_PEDIG_plugin.pedigree_control(pedigree,False,opt.OUTDIR,False) #1st T:1 parent "0" ok.2nd T:relaxed control
    if opt.FIMPUTE: outcome = PHENO_PEDIG_plugin.pedigree_control(pedigree,False,opt.OUTDIR,True) #1st T:1 parent "0" ok.2nd T:strict control
    if not outcome[0]:bomb(outcome[1])
    else:
        if outcome[1]:
            renameped=PHENO_PEDIG_plugin.rename_siredam_pedig(pedigree,opt.OUTDIR)
            if renameped:pedigree=renameped[1]
            warn("Male and female parents in pedigree renamed from UUUUU[...] to 0.\n"+\
                 " "*18+"Using (modified) pedigree file: "+pedigree)
        logit("\n"+" "*14+"[GOOD NEWS]: All pedigree checks ok!")


################################
### OPTION: Plink QC 
################################
if opt.PLINKQC:
    ### Print options in PARAM file
    names = ['Max missings for Individuals','Max missings for SNPs','Minor Allele Frequency',\
           'p(Hardy-Weinberg Equilibrium)','Output file name','Custom (plink) options']
    logit('')
    logit('      ------------------------------------      ')
    logit(' ***  QUALITY CONTROL USING PLINK SOFTWARE  *** ')
    logit('      ------------------------------------     \n')
    
    for i in range(len(names)):
        text = {plinkqc_vals[i]:plinkqc_vals[i],'-9':'NOT REQUIRED'}
        logit('      - %-30s: %-s' % ( names[i],text[plinkqc_vals[i]]))

    ### Do the job
    outcome = PLINK_plugin.runQC_PLINK(PLINK_PATH,pedfile,mapfile,plinkqc_vals,spe)
    if opt.DEBUG:logit('OUTCOME_PLINKQC:\n'+str(outcome)+'\n')
    
    if any('Error:' in string for string in outcome):
        bomb("PLINK QC failed! Please check %s.log to find what went wrong!"% (plinkqc_vals[4]))
    else:
        logit("\n"+" "*14+"[GOOD NEWS]: PLINK QC run OK!. Check         : %s.log " % (plinkqc_vals[4]))
        logit(" "*27+"- Output PLINK (QC'd)           : %s[.map/.ped] " % (plinkqc_vals[4]))
        logit(" "*27+"- PLEASE NOTE: if other analyses are required, the QC'd dataset will be used")

    # Change the ped/map file for the edited version
    pedfile = plinkqc_vals[4]+'.ped'
    mapfile = plinkqc_vals[4]+'.map'

    # If not --save option, delete temp files
    if not opt.SAVEIT:
        logit('\n ***  Deleting TEMP files (add "--save" option if you want to keep all file!)')
        os.system('rm -f '+opt.TMPDIR+'/*')


################################
### OPTION: PEDIGREE checks 
################################
if opt.PEDIG:
    logit('\n      ------------------------------------   ')
    logit(' ***  PEDIGREE CHECKS OF GENOTYPED SAMPLES  ***')
    logit('      ------------------------------------     ')
    logit('\n      '+ptext+'\n')
    ### Do the job
    outcome = PLINK_plugin.runAUTOSOME_PLINK(PLINK_PATH,pedfile,mapfile,opt.TMPDIR,spe)

    if opt.DEBUG:logit('OUTCOME_PLINK_AUTOSOME:\n'+str(outcome)+'\n')
    if any('Error:' in string for string in outcome):
        bomb("PLINK PEDIG - autosome PLINK step - failed! Please check %s.log to find what went wrong!"% (opt.TMPDIR+'/autosomes.log'))

    ### Do the job
    outlist=(opt.OUTDIR+'/PEDIGCHK_pass.txt',opt.OUTDIR+'/PEDIGCHK_fail.txt',opt.OUTDIR+'/PEDIGCHK_bestmatch.txt')
    outcome = PEDIG_plugin.run_pedigree_check(opt.TMPDIR+'/autosomes.ped',pedigree,menderr_thr,skipcouples,outlist,checkall)
    if opt.DEBUG:logit('OUTCOME_PEDIG_CHECK:\n'+str(outcome)+'\n')    

    ### Print the outcome of the analysis
    if outcome[0]:
        names = ['Records read in pedigree file','Genotype records read','   => Used for pedigree check',\
                 '# couples gtyped in pedigree','# couples skipped (from file)','# couples actually controlled',\
                 '          => Mendelcheck PASS','          => Mendelcheck FAIL','     => Total bestmatch found',\
                 ' IN BESTMATCH PROCEDURE ','           Total bestm. checks','         Total SNPs/individual','       Aprox. 1to1 comparisons']
        for i in range(len(names)):logit('      - %-30s: %-s' % (names[i],outcome[i+1]))
        logit("\n"+" "*14+"[GOOD NEWS]: Pedigree checks run OK!")   ## If not stopped, files ok.
        logit(" "*28+"- File containing list of PASS individuals:  %s" % (outlist[0]))
        logit(" "*28+"- File containing list of FAIL individuals:  %s" % (outlist[1]))
        logit(" "*28+"- File containing list of BEST_MATCH individuals:  %s" % (outlist[2]))
    else: bomb(outcome[1])

################################
### OPTION: MEND checks
################################
if opt.MEND:
    logit('\n      -----------------------------------------------------   ')
    logit(' ***  MENDELIAN INHERITANCE CHECKS OF ALL GENOTYPED SAMPLES  ***')
    logit('      -----------------------------------------------------     ')
    ### Do the job
    outcome = PLINK_plugin.runAUTOSOME_PLINK(PLINK_PATH,pedfile,mapfile,opt.TMPDIR,spe)

    if opt.DEBUG:logit('OUTCOME_PLINKMEND_AUTOSOME:\n'+str(outcome)+'\n')
    if any('Error:' in string for string in outcome):
        bomb("PLINK MEND - autosome PLINK step - failed! Please check %s.log to find what went wrong!"% (opt.TMPDIR+'/autosomes.log'))

    ### Do the job
    outcome = MEND_plugin.run_mend_check(opt.TMPDIR+'/autosomes.ped',mend_thr,opt.OUTDIR+'/MENDCHK_lessTHRESHOLD.txt')
    if opt.DEBUG:logit('OUTCOME_MEND_CHECK:\n'+str(outcome)+'\n')

    ### Print the outcome of the analysis
    if outcome[0]:
        logit('      - %-30s: %-s' % ('Total # of genotypes read    ',outcome[3]))
        logit('      - %-30s: %-s' % ('Total # of checks performed  ',outcome[1]))
        logit('      - %-30s: %-s' % ('          => Mendelcheck PASS',outcome[2]))
        logit("\n"+" "*14+"[GOOD NEWS]: Mendelian checks run OK!")   ## If not stopped, files ok.
        logit(" "*28+"- File containing list of individuals below the chosen threshold:  %s" % opt.OUTDIR+'/MENDCHK_lessTHRESHOLD.txt')
    else: bomb(outcome[1])

################################
### OPTION: Plink+R MDS plot
################################
if opt.MDSPLOT:
    logit('      --------------------------------------------     ')
    logit(' ***  MULTI-DIMENSIONAL SCALING PLOT USING PLINK+R  *** ')
    logit('      --------------------------------------------     ')

    ## RUN first part (PLINK)
    logit('\n      Step 1: Calculating MDS coordinates using PLINK')
    outcome = PLINK_plugin.runMDS_PLINK (PLINK_PATH,pedfile,mapfile,opt.TMPDIR,spe)
    if opt.DEBUG:logit('OUTCOME_MDS_PLINK:\n'+str(outcome)+'\n')

    ## CHECK first part (PLINK)
    if any('Error:' in string for string in outcome):
        bomb("PLINK MDS failed! Please check %s.log to find what went wrong!" % (opt.TMPDIR+'/MDSPLOT.log'))
    else:logit("\n"+" "*14+"[GOOD NEWS]: PLINK MDS run OK!")

    ## RUN second part (R)
    if groupop:logit('\n      Step 2: Plotting MDS using R and ggplot2 library (plotting POPULATIONS)')
    else:logit('\n      Step 2: Plotting MDS using R and ggplot2 library (plotting INDIVIDUALS)')
    outcome=PLINK_plugin.runMDS_R (opt.TMPDIR,opt.OUTDIR,groupop,output_name)
    if opt.DEBUG:logit(str(os.system('cat '+opt.TMPDIR+'/mds_plot.Rout')))

    ## CHECK second part (R)
    CHK = open(opt.TMPDIR+'/mds_plot.Rout').readlines()
    if  '[1] "ENDOK"\n' in CHK:
        logit("\n"+" "*14+"[GOOD NEWS]: MDS plot run OK!.")
        if groupop:logit(" "*27+"MDS plot available in         : %s" % outcome)
        else:      logit(" "*27+"MDS plot available in         : %s" % outcome)
    else: bomb("MDS plot failed! Please check "+opt.TMPDIR+"/mds_plot.Rout to find what went wrong!")

    # If not --save option, delete temp files
    if not opt.SAVEIT: 
        logit('\n ***  Deleting TEMP files (add "--save" option if you want to keep all file!)')
        os.system('rm -f '+opt.TMPDIR+'/*')

##########################################
### OPTION: HAPREP (FIRST PART)
##########################################
if opt.HAPREP:
     ### Print options in PARAM file
     logit('\n      --------------------------------------------------   ')
     logit(' ***  PREPARATION OF INPUT FILE FOR SNP2CARRIER SOFTWARE  ***')
     logit('      STEP 1 - Breed chosen: '+opt.HAPREP                     )
     logit('      ----------------------------------------------------   ')
     logit('      -  A) REDUCTION TO CHROMOSOME 19 and IMPUTATION STEP'   )

     ### Read SNPs to keep (By breed)
     if opt.HAPREP=='FLK':alldata=gzip.open(HERE+'/UTILS/ZANARDI_UTILS/haplokeepFLK.txt.gz').readlines()
     elif opt.HAPREP=='BSW':alldata=gzip.open(HERE+'/UTILS/ZANARDI_UTILS/haplokeepBSW.txt.gz').readlines()
     hsnp=dict((alldata[i].strip().split()[0],alldata[i].strip().split()[1:]) for i in range(len(alldata)))
     cromo='19' ### Currently only 1 haplo tested on BTA19

     ### Read map and check consistency
     outcome=HAPLO_plugin.hapchk(mapfile,hsnp,cromo)
     if not outcome[0]:bomb(outcome[1])

     ## First keep only the chromosome we want
     chromo='--chr '+cromo 
     outcome = PLINK_plugin.extract_CHROMO(PLINK_PATH,pedfile,mapfile,opt.TMPDIR,chromo,spe)
     if opt.DEBUG:logit('OUTCOME_PLINK_extract_CHROMO:\n'+str(outcome)+'\n')

     if any('Error:' in string for string in outcome):
         bomb("PLINK extract_CHROMO failed! Please check %s.log to find what went wrong!"% (opt.TMPDIR+'/PLINK_HAPREP'))

     ### Rename ped/map file (all SNPs in chrom 19)
     pedfile=opt.TMPDIR+'/PLINK_HAPREP.ped'
     mapfile=opt.TMPDIR+'/PLINK_HAPREP.map'

     logit('      -  A) REDUCTION TO CHROMOSOME 19 and IMPUTATION STEP - DONE\n')

################################
### OPTION: BEAGLE v.3
################################
if opt.BEAGLE3:
    ### Print options in PARAM file
    names_BG = ['Memory:','Missing Code:','Beagle Output:','Custom (Beagle) options:']
    logit('')
    logit('      ------------------------------      ')
    logit(' ***  BEAGLE v.3 IMPUTATION SOFTWARE  *** ')
    logit('      ------------------------------      ')
    logit('      NOTE: This analysis requires fcGENE and BEAGLE software installed.')

    #CHECK input file. Change Default
    if not beagle_def[0]: warn(' No memory value in PARAMETER file.  Using default value: 2000 Mb')
    if not beagle_def[1]: warn(' No missing value in PARAMETER file. Using default value: "0"')
    
    ## First check consistency of param file
    bgmis=beagle_vals[1]+' '+beagle_vals[1]
    logit('\n      Step 1: Checking presence of user-defined missing values ("'+bgmis+'")')
    found=False
    for i in open(pedfile,'r'):
        string=i.replace('\t',' ').split(' ',6)
        if bgmis in string[6]: found=True
        if found: break
    if not found: warn(' No missing value ("'+bgmis+'") found in ped file - If you expected missings, please:\n'+\
                       ' '*22+'- Set the appropriate missing value in the PARAMETER file (BGMISSING)\n'+\
                       ' '*22+'- Search for the above missing value in the .ped file used ('+pedfile+')\n'+\
                       ' '*22+'- If none of the above works, please contact us.')

    ## RUN first part (fcGENE)
    logit('\n      Step 2: Transforming PLINK data into BEAGLE format (using fcGENE software)')
    outcome = sub.Popen([str(FCGENE_PATH+'./fcgene --ped '+pedfile+' --map '+mapfile+' --oformat beagle'+\
                  ' --out '+opt.TMPDIR+'/PLINK_beagle')],shell=True,stdout=sub.PIPE, stderr=sub.STDOUT).stdout.readlines()
    if opt.DEBUG:logit('OUTCOME_FCGENE(BEAGLE3):\n'+str(outcome)+'\n')

    ## CHECK first part (fcGENE)
    if 'Total time taken for the analysis' in outcome[-1]:logit("\n"+" "*14+"[GOOD NEWS]: fcGENE run OK!")
    else:bomb("fcGENE failed! Please check PLINK_beagle_fcgene.log to find what went wrong!" )

    ## RUN second part (fcGENE)
    logit('\n      Step 3: Imputing (or phasing) data using BEAGLE')
    for i in range(len(names_BG)):
        text = {beagle_vals[i]:beagle_vals[i],'-9':'NOT REQUIRED'}
        logit('      - %-30s %-s' % ( names_BG[i],text[beagle_vals[i]]))
    
    ## RUN Beagle Software
    outcome = sub.Popen([str('java -Xmx'+beagle_vals[0]+'m -jar '+BEAGLE_PATH+\
              'beagle3.jar unphased='+opt.TMPDIR+'/PLINK_beagle.bgl'+' missing='+beagle_vals[1]+' '+\
              beagle_vals[3]+' out='+beagle_vals[2]  )],shell=True,stdout=sub.PIPE, stderr=sub.STDOUT).stdout.readlines()
    if opt.DEBUG:logit('OUTCOME_BEAGLE3:\n'+str(outcome)+'\n')
    
    ## Check Beagle run ok.
    if 'finished' in outcome[-1]:logit("\n"+" "*14+"[GOOD NEWS]: BEAGLE v.3 run OK!")
    else: bomb("BEAGLE failed! Please check "+opt.TMPDIR+"/<beaglenamefile>.log to find what went wrong!")

 ## Conversion Beagle to Plink Format 
    logit('\n      Step 4: Go back to PLINK format')
    outcome = sub.Popen([str(FCGENE_PATH+'./fcgene --bgl '+beagle_vals[2]+\
              '.PLINK_beagle.bgl.phased.gz --oformat plink  --pedinfo '\
                                 +opt.TMPDIR+'/PLINK_beagle_pedinfo.txt --snpinfo '\
              +opt.TMPDIR+'/PLINK_beagle_snpinfo.txt --out '+beagle_vals[2]+'_IMPUTED')],\
               shell=True,stdout=sub.PIPE, stderr=sub.STDOUT).stdout.readlines()
    if opt.DEBUG:logit('OUTCOME_FCGENE_LAST(BEAGLE3):\n'+str(outcome)+'\n')

    ## CHECK Conversion bgl file to Plink first part (fcGENE)
    if 'Total time taken for the analysis' in outcome[-1]:logit("\n"+" "*14+"[GOOD NEWS]: fcGENE Conversion to Plink OK!")
    else:bomb("fcGENE Conversion failed! Please check .log to find what went wrong!" )    

    pedfile=beagle_vals[2]+'_IMPUTED.ped'
    mapfile=beagle_vals[2]+'_IMPUTED.map'

    # If not --save option, delete temp files
    if not opt.SAVEIT: 
        logit('\n ***  Deleting TEMP files (add "--save" option if you want to keep all file!)')
        os.system('rm -f '+opt.TMPDIR+'/*')

################################
### OPTION: BEAGLE v.4 
################################
if opt.BEAGLE4:
    ### Print options in PARAM file
    names_BG = ['Memory:','Missing Code:','Beagle Output:','Custom (Beagle) options:']
    logit('')
    logit('      ------------------------------      ')
    logit(' ***  BEAGLE v.4 IMPUTATION SOFTWARE  *** ')
    logit('      ------------------------------      ')
    logit('      NOTE: This analysis requires PLINK and BEAGLE software installed.')
    
    if not beagle_def[0]: warn(' No Memory value in PARAMETER file.  Using default value: 2000 Mb')
    
    ##check allele type in ped file
    alleletype=0
    lines=open(pedfile,'r').readline()
    if 'B' in lines.replace('\t',' ').strip().split()[6:]: alleletype='A B'
    elif '1' in lines.replace('\t',' ').strip().split()[6:]: alleletype='1 2'

    ##convert allele if necessary 
    if alleletype:
        beagle_ped1=open(opt.TMPDIR+'/BEAGLE4_convert.ped', 'w')
        for line in open(pedfile):
            bre,idx,sire,dame,sex,phe,geno=line.replace('\t',' ').strip().split(' ',6)
            geno_convert=geno.replace("B","C").replace("1","A").replace("2","C")
            beagle_ped1.write('%s %s %s %s %s %s %s \n' % (bre,idx,sire,dame,sex,phe,geno_convert ))
        beagle_ped1.close()
        pedfile=opt.TMPDIR+'/BEAGLE4_convert.ped'        
        warn('PED file contains alleles %s (not ok for Beagle4). The alleles were temporarily converted by Zanardi' % alleletype)

    ##Run a check on map file
    checkid={}
    checkpos={}
    for line in open(mapfile,'r'):
        crom,name,xx,pos=line.replace('\t',' ').strip().split()
        if checkid.has_key(name):
            bomb("Markers with identical ID found ("+name+")! Beagle v.4 would crash, so I'm stopping.\n"+" "*11+\
                 "Please correct (see %s/BEAGLE4_infile.log) and rerun " % (opt.TMPDIR))
        else:checkid[name]=0

    outcome=CHECK_map.snp_position(mapfile,opt.OUTDIR)
    if not outcome[0]:bomb(outcome[1])

    ## Convert PED to VCF using plink 1.9
    logit('\n      Step 1: Convert PED to VCF format using PLINK')
    outcome = PLINK_plugin.vcf_convert (pedfile,mapfile,PLINK_PATH,opt.TMPDIR,spe)
    if opt.DEBUG:logit('OUTCOME_PLINK_VCF_TO_PED(BEAGLE4):\n'+str(outcome)+'\n')

    ## CHECK 
    if any('Error:' in string for string in outcome):
        bomb("VCF convertion failed! Please check PLINK's .log file to find what went wrong!")
    else:
        logit("\n"+" "*14+"[GOOD NEWS]: VCF created OK!")

    ## Do the job
    logit('\n      Step 2: Imputing (or phasing) data using BEAGLE v.4')
    outcome = sub.Popen([str('java -Xmx'+beagle_vals[0]+'m -jar '+BEAGLE_PATH+'beagle4.jar gt='+\
                                 opt.TMPDIR+'/beagle4_infile.vcf '+beagle_vals[3]+' out='+opt.TMPDIR+'/result_beagle4')],\
                            shell=True,stdout=sub.PIPE, stderr=sub.STDOUT).stdout.readlines()
    if opt.DEBUG:logit('OUTCOME_BEAGLE4:\n'+str(outcome)+'\n')

    ## CHECK
    if 'finished' in outcome[-1]:logit("\n"+" "*14+"[GOOD NEWS]: BEAGLE v.4 run OK!")
    else: bomb("BEAGLE v.4 failed! Please check "+opt.TMPDIR+"/result_beagle4.log to find what went wrong!")

    ## Convert VCF to PED using plink1.9
    logit('\n      Step 3: Go back to PLINK format')
    outcome = PLINK_plugin.ped_convert(opt.TMPDIR+'/result_beagle4.vcf.gz',PLINK_PATH,opt.TMPDIR,spe)
    if opt.DEBUG:logit('OUTCOME_PLINK_PED_TO_VCF(BEAGLE4):\n'+str(outcome)+'\n')
    
    ## CHECK
    if any('Error:' in string for string in outcome):
        bomb("PED conversion failed! Please check "+opt.TMPDIR+"/beagle4_imputed.log to find what went wrong!")
    else:logit("\n"+" "*14+"[GOOD NEWS]: PED file created OK!")

    ## Move final files (depending on type of alleles)
    pedfile=beagle_vals[2]+'_IMPUTED.ped'
    mapfile=beagle_vals[2]+'_IMPUTED.map'
    if alleletype:
        logit('\n      Step 4: Converting alleles back to its original status')
        beagle_ped=open(pedfile,'w')
        for line in open(opt.TMPDIR+'/beagle4_imputed.ped'):
            bre,idx,sire,dame,sex,phe,geno=line.strip().split(' ',6)
            if alleletype == 'A B':geno_convert=geno.replace("C","B")
            elif alleletype == '1 2':geno_convert=geno.replace("A","1").replace("C","2")
            beagle_ped.write('%s %s %s %s %s %s %s \n' % (bre,idx,sire,dame,sex,phe,geno_convert ))
        beagle_ped.close()
    else:
        os.system('mv -f '+opt.TMPDIR+'/beagle4_imputed.ped '+pedfile)
    os.system('mv -f '+opt.TMPDIR+'/beagle4_imputed.map '+ mapfile)
    os.system('mv -f '+opt.TMPDIR+'/result_beagle4.log ' + pedfile[:-4]+'.log')

    #Print final log
    logit("\n"+" "*14+"[GOOD NEWS]: BEAGLE v.4 run OK!. Check       : %s " % (pedfile[:-4]+'.log'))
    logit(" "*28+"- Output BEAGLE v.4            : %s[.map/.ped] " % (pedfile[:-4]))

    # If not --save option, delete temp files
    if not opt.SAVEIT: 
        logit('\n ***  Deleting TEMP files (add "--save" option if you want to keep all file!)')
        os.system('rm -f '+opt.TMPDIR+'/*')

#####################################################
### OPTION: FIMPUTE
#####################################################
if opt.FIMPUTE:
    logit('\n      -----------------')
    logit(' ***  FImpute software   ***')
    logit('      -----------------   \n')

    ##check same SNP position in map file ##EZE SPOSTARE QUESTO COME CONTROLLO GENERALE?
    logit('\n      Step 1: Check mapfile')
    outcome=CHECK_map.snp_position(mapfile,opt.OUTDIR)
    if not outcome[0]:bomb(outcome[1])
    chromosome=dict( (map_line.upper().strip().split()[0],0) for map_line in open(mapfile)) 
    if _IMPUTATION_:
        if '0' in chromosome.keys():bomb("Chromosome '0' not allowed in FImpute (required range from 1 to N")

    ##frequency in pedfile ACGT format
    logit('\n      Step 2: Allele frequency calculation using PLINK')
    outcome = FIMPUTE_plugin.allele_freq(PLINK_PATH,pedfile,mapfile,opt.TMPDIR,spe,'freqACGT')
    #check error in plink log
    if any('Error:' in string for string in outcome):
        bomb("PLINK freq failed! Please check "+opt.TMPDIR+"/freqACGT.log to find what went wrong!")
    if opt.DEBUG:logit('Allele Frequency :\n '+opt.TMPDIR+"/freqACGT.frq")

    ## Checking genotypes are not 1/2 (reading first line)
    alleletype=False
    genos=open(pedfile,'r').readline().replace('\t',' ').strip().split()[6:].count('1')
    if genos:alleletype=True

    if not alleletype:
        ## Conversion in 12 format with PLINK
        logit('\n      Step 3: Converting alleles to 1/2 coding format using PLINK')
        outcome = sub.Popen([str(PLINK_PATH+'./plink '+spe+' --recode 12 --ped '+pedfile+' --map '+\
                                     mapfile+' --out '+opt.TMPDIR+'/FIMPUTE_recode12')],\
                                shell=True,stdout=sub.PIPE, stderr=sub.STDOUT).stdout.readlines()
        if any('Error:' in string for string in outcome):
            bomb("PLINK recode12 failed! Please check "+opt.TMPDIR+"/FIMPUTE_recode12.log to find what went wrong!")
        ped_fimpute_recode=opt.TMPDIR+'/FIMPUTE_recode12.ped'
        map_fimpute_recode=opt.TMPDIR+'/FIMPUTE_recode12.map'
    else:
        logit('\n      Step 3: Not converting alleles to 1/2 coding format)')
        ped_fimpute_recode=pedfile
        map_fimpute_recode=mapfile

    if opt.DEBUG:logit('PED_MAP_AFTER_RECODE 1/2:\n Pedfile:\n'+opt.TMPDIR+'/FIMPUTE_recode12.ped'+'\n Mapfile:\n'\
                           +opt.TMPDIR+'/FIMPUTE_recode12.map')

    ## creation file for FImpute
    logit('\n      Step 4: Converting Genotype for Fimpute')
    genotype_out=opt.TMPDIR+'/genotype.FM'
    if output_name:genotype_out=str(opt.TMPDIR+'/genotype'+output_name+'.FM')
    outcome = FIMPUTE_plugin.conversion_PLINK_to_Fimpute(ped_fimpute_recode,genotype_out)
    
    ## creation SNP file for FImpute
    logit('\n      Step 5: Converting SNP map for Fimpute')
    snpinfo_out=opt.TMPDIR+'/snp_info.FM'
    if output_name:snpinfo_out=str(opt.TMPDIR+'/snp_info'+output_name+'.FM')
    outcome = FIMPUTE_plugin.map_convert_FImipute(map_fimpute_recode,snpinfo_out)
    
    ## Process pedigree (if present)
    if fimpute_pedigree:
        logit('\n      Step 6: Pedigree conversion')
        out_pedig=opt.OUTDIR+'/FIMPUTE.pedig'
        if output_name:out_pedig=str(opt.OUTDIR+'/FIMPUTE'+output_name+'.pedig')
        FIMPUTE_plugin.pedig_save(pedigree,out_pedig)
        logit("\n"+" "*14+"[GOOD NEWS]: Pedigree file created OK!")
        logit(" "*27+"- Output Pedigree           : %s " % out_pedig)
        fimpute_othopt=['ped_file = "'+out_pedig+'";\n']
        fimpar[1][1]=fimpar[1][1]+fimpute_othopt[0]
    else:
        logit('\n      ---> Skipping Step 6 (Pedigree not provided)')

    ##creation paramfile for FImpute
    logit('\n      Step 7: Creation Paramater file for FImpute software')
    param_fimpute=opt.TMPDIR+'/param_FImpute.FM'
    if output_name:param_fimpute=str(opt.TMPDIR+'/param_FImpute'+output_name+'.FM')
    output_folder=opt.OUTDIR+'/output_FImpute'
    if output_name:output_folder=str(opt.OUTDIR+'/OUTPUT_FImpute'+output_name)
    outcome = FIMPUTE_plugin.param_FImpute(genotype_out,snpinfo_out,opt.TMPDIR,output_folder,param_fimpute,fimpar)

    ## Do the job
    logit('\n      Step 8: Imputing data using FImpute v.4')
    outcome = sub.Popen([str(FIMPUTE_PATH+'./FImpute '+param_fimpute+' -o')],\
                           shell=True,stdout=sub.PIPE, stderr=sub.STDOUT).stdout.readlines()

    ## Check FImpute run ok.
    if 'Completed' in outcome[-1]:logit("\n"+" "*14+"[GOOD NEWS]: FImpute run OK!")
    else:
        output_error=open(opt.TMPDIR+'/FImpute_error.log','w')
        output_error.write('%s'%(' '.join(outcome)))
        bomb("FImpute failed! Please check "+opt.TMPDIR+"/FImpute_error.log to find what went wrong!")
    if opt.DEBUG:logit('OUTCOME_FIMPUTE:\n'+str(outcome)+'\n')

    ##rename map file
    final_map=opt.OUTDIR+'/FIMPUTE.map'
    if output_name:final_map=str(opt.OUTDIR+'/FIMPUTE'+output_name+'.map')

    #converted FImpute file to plink format
    logit('\n      Step 9: Convert FImput output in PLINK format (.ped/.map)')
    final_ped=opt.OUTDIR+'/FIMPUTE.ped'
    if output_name:final_ped=str(opt.OUTDIR+'/FIMPUTE'+output_name+'.ped')
    outcome = FIMPUTE_plugin.conversion_Fimpute_to_PLINK(output_folder+'/genotypes_imp.txt',opt.TMPDIR,final_ped,pedfile,output_folder,final_map)

    if outcome:
        logit("\n"+" "*14+"[GOOD NEWS]: Imputed Ped file created OK!")
        logit(" "*27+"- Output FImpute [.ped/.map]  : %s " % final_ped.replace('.ped','[.ped/.map]'))
        logit(" "*27+"- Output FImpute (folder)     : %s " % output_folder)
        ##change new .ped/.map file
        pedfile=final_ped
        mapfile=final_map

    # If not --save option, delete temp files
    if not opt.SAVEIT:
        os.system('rm -f '+opt.TMPDIR+'/*')

######################################
### OPTION: HAPLOTYPE PREP (PART TWO) 
######################################
if opt.HAPREP:
    logit('\n      --------------------------------------------------   ')
    logit(' ***  PREPARATION OF INPUT FILE FOR SNP2CARRIER SOFTWARE  ***')
    logit('      STEP 2 - Breed chosen: '+opt.HAPREP                     )
    logit('      ----------------------------------------------------   ')
    logit('      -  B) REDUCTION TO SELECTED SNPs and OUTPUT PREPARATION')
    
    ### Currently writing a list of file as --snps function is currently not working OK in PLINK1.9 (June 2015)
    retain=open(opt.TMPDIR+'/HAPLO_small.txt','w')
    retain.write('\n'.join(hsnp.keys())+'\n');retain.close()
    chromo='--d : --extract '+opt.TMPDIR+'/HAPLO_small.txt'
    
    ### Keep the 18 SNPs of interest only
    outcome = PLINK_plugin.extract_CHROMO(PLINK_PATH,pedfile,mapfile,opt.TMPDIR,chromo,spe)

    ### Rename ped/map to (short version of 18 SNPs) IMPUTED genotypes
    pedfile=opt.TMPDIR+'/PLINK_HAPREP.ped'
    mapfile=opt.TMPDIR+'/PLINK_HAPREP.map'

    ### Check outcome
    if opt.DEBUG:logit('OUTCOME_PLINK_extract_CHROMO2:\n'+str(outcome)+'\n')
    if any('Error:' in string for string in outcome):
        bomb("PLINK extract_CHROMO failed! Please check %s.log to find what went wrong!"% (opt.TMPDIR+'/PLINK_HAPREP'))

    ### Create input file for HAPLOTYPE testing
    outcome=HAPLO_plugin.hapconvert(opt.TMPDIR,opt.OUTDIR,output_name,mapfile,pedfile,hsnp,opt.HAPREP)
    if outcome[0]:
        logit("\n"+" "*14+"[GOOD NEWS]:\n")
        logit(' '*14+' - SNP2CARRIER input file created correctly!! YAAAAYYYY!!!')
        logit(' '*14+' - Please visit =>  https://stebif68.shinyapps.io/EzeApp  to run your analysis.')
        logit(' '*14+' - Use: '+opt.OUTDIR+'/HAPREP'+str(output_name)+'.txt as input file\n'+' '*14+' - NOW EXITING ZANARDI...')
    else:bomb(outcome[1])

    # If not --save option, delete temp files
    if not opt.SAVEIT:
        logit('\n ***  Deleting TEMP files (add "--save" option if you want to keep all file!)')
        os.system('rm -f '+opt.TMPDIR+'/*')
    sys.exit()

################################
### OPTION: ROH & FROH
################################
if opt.ROH or opt.FROH:
    logit('')
    logit('       --------------------')
    logit('  ***  RUNS OF HOMOZYGOSITY  *** ')
    logit('       --------------------\n')

    ## Conversion in 12 format with PLINK
    logit('\n      Step 1: Converting allele format to 1/2 (required by the pgm) using PLINK')

    outcome = sub.Popen([str(PLINK_PATH+'./plink '+spe+' --recode 12 --ped '+pedfile+' --map '+\
                                 mapfile+' --out '+opt.TMPDIR+'/ROH_recode12')],\
                            shell=True,stdout=sub.PIPE, stderr=sub.STDOUT).stdout.readlines()

    #check error in plink log
    if any('Error:' in string for string in outcome):
        bomb("PLINK recode12 failed! Please check Zanardi.log to find what went wrong!")
    if opt.DEBUG:logit('PED_MAP_AFTER_RECODE 1/2:\n Pedfile:\n'+opt.TMPDIR+'/ROH_recode12.ped'+'\n Mapfile:\n'\
                           +opt.TMPDIR+'/ROH_recode12.map')
    
    ## Read mapfile
    logit('\n      Step 2: Importing Mapfile')
    outcome_map = ROH_plugin.read_map(opt.TMPDIR+'/ROH_recode12.map',nchrom)

    ## Search Runs of Homozygosity
    logit('\n      Step 3: Search Runs Of Homozygosity ')
    outcome = ROH_plugin.ROH(opt.TMPDIR+'/ROH_recode12.ped',outcome_map,roh_vals)

    logit('\n ==>  %-30s' % 'ROH Parameters:')
    logit('      - %-20s %-s' % ( 'Minimum SNP in ROH:  ',roh_vals[0]))
    logit('      - %-20s %-s' % ( 'Minimum ROH length:  ',roh_vals[3]))
    logit('      - %-20s %-s' % ( 'Missings allowed:    ',roh_vals[1]))
    logit('      - %-20s %-s' % ( 'Heterozygous allowed:',roh_vals[2]))
    logit('      - %-20s %-s' % ( 'Name file output:    ',roh_vals[4]))

    n_breed=outcome[0]
    if  outcome[1]:bomb("ROH search failed! Please check ROH.log file to find what went wrong!")
    else:logit("\n"+" "*14+"[GOOD NEWS]: ROH file created OK!")
    
    if opt.DEBUG:logit('ROH_FILE: '+roh_vals[4])    
    
    ## Inbreeding coefficient
    if opt.FROH:
        logit('\n      Step 4: Calculating Inbreeding Coefficients ')
        if output_name:inbreeding_output=opt.OUTDIR+'/ROH_inbreeding'+str(output_name)+'.txt'
        else:inbreeding_output=opt.OUTDIR+'/ROH_inbreeding.txt'
        outcome = ROH_plugin.inbreeding_ROH(roh_vals[4],\
                          opt.TMPDIR+'/ROH_recode12.map',inbreeding_output,nchrom )    
        if inbreeding_output:logit("\n"+" "*14+"[GOOD NEWS]: File created OK!")
        else:bomb("Creation failed!")
    else:
        logit('\n      ---> Skipping Step 4 (Inbreeding Coefficient)')

    ## PLOT Runs of Homozygosity
    if opt.ROH:
        logit('\n      Step 5: Counting SNPs within ROHs ')
        outcome = ROH_plugin.snp_inside_ROH(roh_vals[4],opt.TMPDIR+'/ROH_recode12.map',opt.TMPDIR+'/SNP_inside_roh.txt',n_breed,nchrom)

        if not opt.TMPDIR+'/SNP_inside_roh.txt':bomb("Creation failed!")
        
        logit('\n      Step 6: Creating ROH plot ')
        #roh_vals[5]=roh_vals[5].replace('.txt','').replace('.csv','')

        outcome=ROH_plugin.plot_ROH(opt.TMPDIR,opt.OUTDIR,output_name)
        if opt.DEBUG:logit(str(os.system('cat '+opt.TMPDIR+'/roh_plot.Rout')))
        
        CHK = open(opt.TMPDIR+'/roh_plot.Rout').readlines()
        if  '[1] "ENDOK"\n' in CHK:
            logit("\n"+" "*14+"[GOOD NEWS]: ROH plot run OK!")
            logit(" "*21+"ROH plot available in: %s%s" % (opt.OUTDIR,outcome))
        else: bomb("ROH plot failed! Please check "+opt.TMPDIR+"/roh_plot.Rout to find what went wrong!")
    

#####################################################
### OPTION: ADMIXTURE
#####################################################
if opt.ADMIXTURE:
    logit('\n      ------------------------------------------------------      ')
    logit(' ***  POPULATION STRUCTURE ANALYSIS USING ADMIXTURE SOFTWARE  ***')
    logit('      ------------------------------------------------------    \n')
    logit('\n ==>  %-30s' % 'Admixture Parameters:')
    logit('      - %-20s %-s' % ( 'K Value :    ',admixture_vals[0]))
    logit('      - %-20s %-s' % ( 'Number of threads : ',admixture_vals[1]))
    logit('      - %-20s %-s' % ( 'Cross-Validation :  ',admixture_vals[2]))

    ## Conversion in 12 format with PLINK
    logit('\n      Step 1: Converting allele format to 1/2 (required by the pgm) using PLINK')
    outcome = sub.Popen([str(PLINK_PATH+'./plink '+spe+' --recode 12 --ped '+pedfile+' --map '+\
                                 mapfile+' --out '+opt.TMPDIR+'/Admixture'+output_name+"_K")],\
                            shell=True,stdout=sub.PIPE, stderr=sub.STDOUT).stdout.readlines()

    #check error in plink log
    if any('Error:' in string for string in outcome):
        bomb("PLINK recode12 failed! Please check Zanardi.log to find what went wrong!")
    if opt.DEBUG:logit('PED_MAP_AFTER_RECODE 1/2:\n Pedfile:\n'+opt.TMPDIR+'/Admixture'+output_name+'_K.ped'+\
                           '\n Mapfile:\n'+opt.TMPDIR+'/Admixture'+output_name+'_K.map')
    
    ## Read mapfile
    logit('\n      Step 2: Admixture Analysis')
    os.system('rm -f '+opt.TMPDIR+'/Admixture_CVvalues'+output_name+'.txt')
    CV_admixture_output=open(opt.TMPDIR+"/Admixture_CVvalues"+output_name+".txt",'a')

    for k in range(2,int(admixture_vals[0])+1):
        outcome = sub.Popen([str(ADMIXTURE_PATH+'admixture  '+opt.TMPDIR+'/Admixture'+output_name+'_K.ped '+str(k)+\
                                     ' --cv='+admixture_vals[2]+' -j'+admixture_vals[1] )], \
                                shell=True,stdout=sub.PIPE, stderr=sub.STDOUT).stdout.readlines()
                 
        #Check Error
        if not any('Writing' in string for string in outcome): 
            admixture_output=open(opt.TMPDIR+"/admixture"+output_name+".log",'w')
            for line in outcome:admixture_output.write(line)
            admixture_output.close()
            bomb("ADMIXTURE analysis failed! Please check "+opt.TMPDIR+"/admixture"+output_name+".log to find what went wrong!")
        else:
            os.system('mv -f *.P *.Q '+opt.OUTDIR)
            cv = ([ line for line in outcome  if 'CV' in line.strip().split() ])
            CV_admixture_output.write((str(cv[0])))
    CV_admixture_output.close()

    ##PLOT CV admixture
    logit('\n      Step 3: CV Admixture Plot')
    outcome = ADMIXTURE_plugin.CV_ADM_plot(opt.TMPDIR,opt.OUTDIR,'/Admixture_CVvalues'+output_name+'.txt',output_name)
    CHK = open(opt.TMPDIR+'/CV_admixture_plot.Rout').readlines()
    if  '[1] "ENDOK"\n' in CHK:
        logit(" "*17+"CV Admixture plot available in : %s" % opt.OUTDIR+'/Admixture_CVplot'+output_name+'.pdf' )
    else: bomb("CV Admixture plot failed! Please check "+opt.TMPDIR+"/CV_admixture_plot.Rout to find what went wrong!")    

    ##PLOT Admixture
    logit('\n      Step 4: Admixture Plot')
    os.system("cut -b 1-300  "+opt.TMPDIR+'/Admixture'+output_name+"_K.ped  | awk '{print $1,$2}' > "+opt.TMPDIR+'/Admixture_names.txt')

    outcome = ADMIXTURE_plugin.admixture_plot(opt.TMPDIR,opt.OUTDIR,output_name,int(admixture_vals[0]))

    CHK = open(opt.TMPDIR+'/Admixture_plot.Rout').readlines()
    if  '[1] "ENDOK"\n' in CHK:
        logit(" "*17+"Admixture plot available in    : %s" % opt.OUTDIR+'/Admixture_BARplot'+output_name+'.pdf')
    else: bomb("Admixture plot failed! Please check "+opt.TMPDIR+"/Admixture_plot.Rout to find what went wrong!")    

    # If not --save option, delete temp files
    if not opt.SAVEIT: 
        logit('\n ***  Deleting TEMP files (add "--save" option if you want to keep all file!)')
        os.system('rm -f '+opt.TMPDIR+'/*')


if not opt.SAVEIT: 
    logit('\n ***  Deleting TEMP folder (add "--save" option if you want to keep this folder!)')
    os.system('rm -rf '+opt.TMPDIR)
else:logit('\n ***  Keeping TEMP folder as requested ("--save" option)')
logit(' ***  Final files available in: '+opt.OUTDIR)

if opt.VER:print "\n\nBAZINGA! I've done my stuff."
logit('*'*81+'\n ==> PROGRAM ENDED SUCCESSFULLY: '+time.strftime("%B %d, %Y - %l:%M%p %Z")+' <==\n'+'*'*81)
