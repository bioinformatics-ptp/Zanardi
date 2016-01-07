import subprocess as sub

############################
#### PLINK MERGE plugin ####
############################
def runMERGE_PLINK(PLINK_PATH,pedfile,mapfile,TMPDIR,species,OUTDIR):
    temp_merge=open(TMPDIR+'/merge_list.txt','w')
    for num in range(len(pedfile[1:])):temp_merge.write('%s %s\n' % (pedfile[num+1],mapfile[num+1]))
    temp_merge.close()
    std = sub.Popen([str(PLINK_PATH+'./plink --noweb '+species+' --ped '+pedfile[0]+' --map '+mapfile[0]+\
          ' --merge-list '+TMPDIR+'/merge_list.txt --recode --out '+OUTDIR+'/PLINK_MERGED')],shell=True,\
          stdout=sub.PIPE, stderr=sub.STDOUT).stdout.readlines()
    return std

##########################
#### PLINK STD plugin #### (for when sep is not blank space)
##########################
def runSTD_PLINK(PLINK_PATH,pedfile,mapfile,TMPDIR,species,OUTDIR):
    std = sub.Popen([str(PLINK_PATH+'./plink --noweb '+species+' --ped '+pedfile+' --map '+mapfile+\
          ' --recode --out '+OUTDIR+'/PLINK_STND')],shell=True,\
          stdout=sub.PIPE, stderr=sub.STDOUT).stdout.readlines()
    return std

#########################
#### PLINK QC plugin ####
#########################
def runQC_PLINK(PLINK_PATH,pedfile,mapfile,plinkqc_vals,species):
    plinkopts=[]
    qcopts=[' --mind',' --geno',' --maf',' --hwe',' --out','']
    for i in range(len(plinkqc_vals)):
        if plinkqc_vals[i]=='-9':continue #Skip params not provided
        plinkopts.append(qcopts[i]+' '+plinkqc_vals[i]) 
    std = sub.Popen([str(PLINK_PATH+'./plink --noweb '+species+' --ped '+pedfile+' --map '+mapfile+\
          ' '.join(plinkopts)+' --recode')],shell=True,stdout=sub.PIPE, stderr=sub.STDOUT).stdout.readlines()
    return std

#########################
#### PLINK AUTOSOMES ####
#########################
def runAUTOSOME_PLINK(PLINK_PATH,pedfile,mapfile,TMPDIR,species):
    std = sub.Popen([str(PLINK_PATH+'./plink --noweb '+species+' --ped '+pedfile+' --map '+mapfile+\
          ' --autosome --recode --out '+TMPDIR+'/autosomes')],shell=True,stdout=sub.PIPE, stderr=sub.STDOUT).stdout.readlines()
    return std

##############################
#### MDS PLINK plugin (1) ####
##############################
def runMDS_PLINK(PLINK_PATH,pedfile,mapfile,TMPDIR,species):
    std = sub.Popen([str(PLINK_PATH+'./plink --noweb '+species+' --ped '+pedfile+' --map '+mapfile+\
          ' --autosome --cluster --mds-plot 2 --out '+TMPDIR+'/MDSPLOT')],shell=True,stdout=sub.PIPE, stderr=sub.STDOUT)\
          .stdout.readlines()
    return std

##########################
#### MDS R plugin (2) ####
##########################
def runMDS_R(TMPDIR,OUTDIR,groupop,fileoutput):
    OUTmds=TMPDIR+"/MDSPLOT"
    mds=open(TMPDIR+'/mds_plot.R','w')
    mds.write("setwd('"+TMPDIR+"')\n")
    mds.write("library('ggplot2')\n")
    mds.write("data<-read.table('"+OUTmds+".mds',header=T)\n")
    mds.write("names(data)<-c('Population','ID','sol','PC1','PC2')\n")
    if not groupop:
        pdf_file_name=OUTDIR+'/MDS_plot'+fileoutput+'.pdf'
        mds.write("pdf('"+pdf_file_name+"')\n")
        mds.write("ggplot(data,aes(x=PC2,y=PC1,color=Population,title='MDS plot - Individuals by Population')) " +\
                  "+ geom_point(size=3) + guides(col=guide_legend(ncol=ceiling(length(unique(data$Population))/15)))\n")
    else:
        pdf_file_name=OUTDIR+'/MDS_PLOT_Pops'+fileoutput+'.pdf'
        mds.write("pdf('"+pdf_file_name+"')\n")
        mds.write("POPS<-as.data.frame(tapply(data$PC1,data$Population,mean))\n")
        mds.write("POPS$PC2<-tapply(data$PC2,data$Population,mean)\n")
        mds.write("POPS$Pop<-row.names(POPS)\n")
        mds.write("names(POPS)<-c('PC1','PC2','Population')\n")
        mds.write("ggplot(POPS,aes(x=PC2,y=PC1,color=Population,title='MDS plot - Avg. by population')) "+ \
                  "+ geom_point(size=3) + guides(col=guide_legend(ncol=ceiling(length(unique(POPS$Population))/15)))\n")
    mds.write('dev.off()\n')
    mds.write("print('ENDOK')\n")
    mds.close()

    ### LAUNCH PGM
    sub.call(["R CMD BATCH %s/mds_plot.R %s/mds_plot.Rout " % (TMPDIR,TMPDIR)],shell=True)

    return (pdf_file_name)

######################
#### PLINK to VCF ####
######################
def vcf_convert(pedfile,mapfile,PLINK_PATH,TMPDIR,species):
    std = sub.Popen([str(PLINK_PATH+'./plink --noweb '+species+' --ped '+pedfile+' --map '+mapfile+\
          ' --recode vcf --out '+TMPDIR+'/beagle4_infile')],shell=True,stdout=sub.PIPE, stderr=sub.STDOUT)\
          .stdout.readlines()
    return std

#######################
#### VCF to PLINK  ####
#######################
def ped_convert(vcffile,PLINK_PATH,TMPDIR,species):
    std = sub.Popen([str(PLINK_PATH+'./plink --noweb '+species+' --vcf '+vcffile+' --recode --out '+TMPDIR+'/beagle4_imputed')],\
                        shell=True,stdout=sub.PIPE, stderr=sub.STDOUT).stdout.readlines()
    return std


############################
#### PLINK extractCROMO ####
############################
def extract_CHROMO(PLINK_PATH,pedfile,mapfile,TMPDIR,chromo,species):
    std = sub.Popen([str(PLINK_PATH+'./plink --noweb '+species+' --ped '+pedfile+' --map '+mapfile+\
          ' '+chromo+' --recode --out '+TMPDIR+'/PLINK_HAPREP')],shell=True,stdout=sub.PIPE, stderr=sub.STDOUT).stdout.readlines()
    return std
