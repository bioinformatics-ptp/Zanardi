import subprocess as sub

def runit(usropt,UTILDIR,SyS,HERE,DEBUG):
    ### List software to download - check if option is within options

    DOWN_opts = ['plink','fcgene','beagle3','beagle4','admixture','fimpute']
    for prog in usropt:
        if not prog in DOWN_opts:return(False,"Invalid software request ("+prog+") -> See Zanardi's manual for accepted values")

    ### Build the command launching download.app
    command=UTILDIR+"/ZANARDI_UTILS/./download_app.sh "+SyS+" "+ UTILDIR +' '+ HERE+' '+\
            ' '.join(usropt)+' &> ' + UTILDIR + "/DOWNLOAD.log"
    if DEBUG:print('COMMAND LAUNCHED:'+str(command))
    sub.Popen([command],shell = True, stdout = sub.PIPE,stderr = sub.STDOUT, executable='/bin/bash').stdout.readline().strip()

    ### check download.log file and build returns
    CHK=sub.Popen(["tail -1 "+UTILDIR+"/DOWNLOAD.log"], shell=True,stdout=sub.PIPE, \
        stderr=sub.STDOUT).stdout.readline().strip()
    if 'BAZINGA' in CHK:
        return(True,"\n"+" "*14+"[GOOD NEWS]: Software downloaded OK! Your PARAMETER.txt file has been updated!"+\
                    "\n"+" "*29+"- The required software is available in "+UTILDIR+\
                    "\n"+" "*29+"- Please check download log: "+UTILDIR+"/DOWNLOAD.log"+\
                    "\n"+" "*29+"  (now you should re-run Zanardi with the required analysis!")
    else: return(False,"Download failed! Please check "+UTILDIR+"/DOWNLOAD.log to find what went wrong!")
