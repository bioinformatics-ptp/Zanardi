#!/bin/bash
set -e
###############################################################################
### Simple bash script to download binary (or source codes) of the software ###
### required to run Zanardi options.                                        ###
### History: Feb 2015 - Ezequiel L. Nicolazzi: Original version             ###
###                                                                         ###
### For bug report/comments: ezequiel.nicolazzi@tecnoparco.org              ###
#######################################################ELN#####################

### Variables passed by Zanardi
SYS=$1       # This is the OS
WHERE=$2     # This is the path where to download the files
SOURCE=$3    # This is the path where Zanardi and its param file are
DWNL=${@:3}  # This is the list of software to download

####################################################################################
###                            List of download links                            ###
####################################################################################
if [ $SYS == 'Linux' ]; then
    method=wget
    opt=-O
    PLINK=http://www.cog-genomics.org/static/bin/plink170724/plink_linux_x86_64.zip
    FCGENE=http://sourceforge.net/projects/fcgene/files/latest/download
    ADMIXTURE=https://www.genetics.ucla.edu/software/admixture/binaries/admixture_linux-1.23.tar.gz
    FIMPUTE=http://www.aps.uoguelph.ca/~msargol/fimpute/FImpute_Linux.zip
else
    method=curl
    opt=-O
    PLINK=http://www.cog-genomics.org/static/bin/plink170724/plink_mac.zip
    FCGENE=http://heanet.dl.sourceforge.net/project/fcgene/fcgene-1.0.7.tar.gz
    ADMIXTURE=https://www.genetics.ucla.edu/software/admixture/binaries/admixture_macosx-1.3.0.tar.gz
    FIMPUTE=http://www.aps.uoguelph.ca/~msargol/fimpute/FImpute_Mac.zip
fi

BEAGLE3=https://faculty.washington.edu/browning/beagle/beagle.jar
BEAGLE4=https://faculty.washington.edu/browning/beagle/beagle.r1399.jar

####################################################################################

#### Routine for mac system (own programs)
for pgm in $DWNL  
do
    if [ $pgm == 'plink' ]; then
	pgm_folder=$WHERE/PLINK 
	echo;echo " ==> Downloading PLINK in $WHERE/PLINK from $PLINK using $method";echo
	if [ $SYS == 'Linux' ];then
	    mkdir -p $pgm_folder && $method $PLINK $opt $pgm_folder/plink.zip ||  exit
	else
	    mkdir -p $pgm_folder && cd $pgm_folder && $method $opt $PLINK && mv $pgm_folder/plink_mac.zip $pgm_folder/plink.zip ||  exit
	fi
	echo;echo " ==> Uncompressing the file .zip";echo
	cd $pgm_folder && unzip -o plink.zip || exit && dirplink=$pgm_folder         #plink1.9
	#chmod 755 $dirplink/plink    #plink1.07
	chmod -R 755 $dirplink        #plink1.9
	changepath=`grep PGM_PLINK $SOURCE/PARAMFILE.txt`
	sed -i -e "s#$changepath#PGM_PLINK=$dirplink#" $SOURCE/PARAMFILE.txt
	echo;echo " ==> PLINK downloaded and uncompressed successfully, path updated in PARAMETER.txt" 
    elif [ $pgm == 'fcgene' ]; then
	pgm_folder=$WHERE/FCGENE
	echo;echo " ==> Downloading fcGENE in $WHERE/FCGENE from $FCGENE using $method";echo
        if [ $SYS == 'Linux' ];then
 	    mkdir -p $pgm_folder && $method $FCGENE $opt $pgm_folder/fcgene.tar.gz || exit
	else
            mkdir -p $pgm_folder && cd $pgm_folder && $method $opt $FCGENE && mv $pgm_folder/fcgene-1.0.7.tar.gz $pgm_folder/fcgene.tar.gz ||  exit
	fi
	echo;echo " ==> Uncompressing the file .tar.gz";echo
	cd $pgm_folder && tar -zxvf fcgene.tar.gz  && dirfcgene=$pgm_folder/`ls -d */`
	chmod 755 $dirfcgene/fcgene
	changepath=`grep PGM_FCGENE $SOURCE/PARAMFILE.txt`
	sed -i -e "s#$changepath#PGM_FCGENE=$dirfcgene#" $SOURCE/PARAMFILE.txt
	echo;echo " ==> FCGENE downloaded and uncompressed successfully, path updated in PARAMETER.txt" 
    elif [ $pgm == 'beagle3' ]; then
	pgm_folder=$WHERE/BEAGLE
	echo;echo " ==> Downloading BEAGLE v.3 in $WHERE/BEAGLE from $BEAGLE using $method";echo
        if [ $SYS == 'Linux' ];then
	    mkdir -p $pgm_folder && $method $BEAGLE3 $opt $pgm_folder/beagle3.jar || exit
	else
            mkdir -p $pgm_folder && cd $pgm_folder && $method $opt $BEAGLE3 && mv $pgm_folder/beagle.jar $pgm_folder/beagle3.jar ||  exit
	fi	    
	chmod 755 $pgm_folder/beagle3.jar
	changepath=`grep PGM_BEAGLE3 $SOURCE/PARAMFILE.txt`
	sed -i -e "s#$changepath#PGM_BEAGLE3=$pgm_folder#" $SOURCE/PARAMFILE.txt
	echo;echo " ==> BEAGLE v.3 downloaded successfully, path updated in PARAMETER.txt" 
    elif [ $pgm == 'beagle4' ]; then
	pgm_folder=$WHERE/BEAGLE
	echo;echo " ==> Downloading BEAGLE v.4 in $WHERE/BEAGLE from $BEAGLE using $method";echo
	if [ $SYS == 'Linux' ];then
	    mkdir -p $pgm_folder && $method $BEAGLE4 $opt $pgm_folder/beagle4.jar || exit
	else
	    mkdir -p $pgm_folder && cd $pgm_folder && $method $opt $BEAGLE4 && mv $pgm_folder/beagle.r1399.jar $pgm_folder/beagle4.jar ||  exit
	fi
	chmod 755 $pgm_folder/beagle4.jar
	changepath=`grep PGM_BEAGLE4 $SOURCE/PARAMFILE.txt`
	sed -i -e "s#$changepath#PGM_BEAGLE4=$pgm_folder#" $SOURCE/PARAMFILE.txt
	echo;echo " ==> BEAGLE v.4 downloaded successfully, path updated in PARAMETER.txt" 
    elif [ $pgm == 'admixture' ]; then
	pgm_folder=$WHERE/ADMIXTURE
	echo;echo " ==> Downloading ADMIXTURE in $WHERE/ADMIXTURE from $ADMIXTURE using $method";echo
        if [ $SYS == 'Linux' ];then
	    mkdir -p $pgm_folder && $method $ADMIXTURE $opt $pgm_folder/admixture.tar || exit
        else
            mkdir -p $pgm_folder && cd $pgm_folder && $method $opt $ADMIXTURE && mv $pgm_folder/admixture_macosx-1.3.0.tar.gz $pgm_folder/admixture.tar ||  exit
        fi
	echo;echo " ==> Uncompressing the file .tar";echo
	cd $pgm_folder && tar -xvf admixture.tar  && diradmixture=$pgm_folder/`ls -d */`
	chmod -R 755 $diradmixture
	changepath=`grep PGM_ADMIXTURE $SOURCE/PARAMFILE.txt`
	sed -i -e "s#$changepath#PGM_ADMIXTURE=$diradmixture#" $SOURCE/PARAMFILE.txt
	echo;echo " ==> ADMIXTURE downloaded successfully, path updated in PARAMETER.txt" 
    elif [ $pgm == 'fimpute' ]; then
	pgm_folder=$WHERE/FIMPUTE
	echo;echo " ==> Downloading FIMPUTE in $WHERE/FIMPUTE from $FIMPUTE using $method";echo
	if [ $SYS == 'Linux' ];then
	    mkdir -p $pgm_folder && $method $FIMPUTE $opt $pgm_folder/fimpute.zip ||  exit
	else
            mkdir -p $pgm_folder && cd $pgm_folder && $method $opt $FIMPUTE && mv $pgm_folder/FImpute_Mac.zip $pgm_folder/fimpute.zip ||  exit
	fi
	echo;echo " ==> Uncompressing the file .zip";echo

	if [ $SYS == 'Linux' ];then
	    cd $pgm_folder && unzip -o fimpute.zip || exit 
	    mv $pgm_folder/FImpute_Linux $pgm_folder/FImpute
	else
	    cd $pgm_folder && unzip -o fimpute.zip || exit 
	    mv $pgm_folder/FImpute_Mac $pgm_folder/FImpute
	    mv $pgm_folder/FImpute/FImpute_mac $pgm_folder/FImpute/FImpute
	fi
	dirfimpute=$pgm_folder/FImpute
	chmod -R 755 $dirfimpute        
	changepath=`grep PGM_FIMPUTE $SOURCE/PARAMFILE.txt`
	sed -i -e "s#$changepath#PGM_FIMPUTE=$dirfimpute#" $SOURCE/PARAMFILE.txt
	echo;echo " ==> FImpute downloaded and uncompressed successfully, path updated in PARAMETER.txt" 
    fi
done
rm -f $SOURCE/PARAMFILE.txt-e
echo "BAZINGA"
