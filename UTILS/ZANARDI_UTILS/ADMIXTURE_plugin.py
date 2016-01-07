import subprocess as sub

#plot cross-validation
def CV_ADM_plot(TMPDIR,OUTDIR,data_file,name_file):
    if name_file:fileoutput='_'+name_file
    else:fileoutput=name_file
    OUT_adm=TMPDIR+str(data_file)
    adm=open(TMPDIR+'/CV_admixture_plot.R','w')
    adm.write("setwd('"+TMPDIR+"')\n")
    adm.write("library('ggplot2')\n")
    adm.write("data_admixture<-read.table('"+OUT_adm+"',header=F)\n")
    adm.write("pdf('"+OUTDIR+"/Admixture_CVplot"+name_file+".pdf',height=12, width=20)\n")
    adm.write("data_admixture$V5=as.integer(gsub('"'\\\(\\\K='"','',gsub('.{2}$','',data_admixture$V3))) \n")
    adm.write(" grafico=ggplot(data = data_admixture, aes(x=V5,y=V4) ) + geom_line() + geom_point() + xlab('"'K '"') + ylab('"'Cross-validation error'"') \n")
    adm.write("print(grafico)  \n dev.off() \n")
    adm.write("print('ENDOK')\n")
    adm.close()

    ##LAUNCH PGM
    sub.call(["R CMD BATCH %s/CV_admixture_plot.R %s/CV_admixture_plot.Rout " % (TMPDIR,TMPDIR)],shell=True)

    return (fileoutput)

# BAR plot
def admixture_plot(TMPDIR,OUTDIR,fileoutput,k_value):
    adm=open(TMPDIR+'/Admixture_plot.R','w')
    adm.write("library('ggplot2')\n")
    adm.write("pdf('"+OUTDIR+"/Admixture_BARplot"+fileoutput+".pdf',height=12, width=20) \n")
    adm.write("for (k in seq(2,"+str(k_value)+")){ \n")
    adm.write("data=read.table(paste('"+OUTDIR+"/Admixture"+fileoutput+"_K.',k,'.Q',sep=''),stringsAsFactors = TRUE)\n ")
    adm.write("names=read.table('"+TMPDIR+"/Admixture_names.txt')\n")
    adm.write("names$V3<-paste(names$V1,names$V2,sep='_')\n")
    adm.write("barplot(t(as.matrix(data)), col=rainbow(k), main=paste('K = ',k,sep=' '),xlab='"'Individual #'"', ylab='"'Ancestry'"', border=NA,names.arg=names$V3, cex.names=.5, las=2) \n")
    adm.write("}\n dev.off()\n")
    adm.write("print('ENDOK')\n")
    adm.close()

    ##LAUNCH PGM
    sub.call(["R CMD BATCH %s/Admixture_plot.R %s/Admixture_plot.Rout " % (TMPDIR,TMPDIR)],shell=True)

    return (fileoutput)


