#!/bin/bash
. /seq/vgb/jason/scripts/bashFunctions.sh

#snpeffDir="/seq/vgb/software/snpEff_4.2/"
snpeffDir="/seq/vgb/swofford/software/snpEff/"
description="runs snpeff using gridengine on the specified vcf file"
argumentDef=()
argumentDef+=("-vcf" vcfFile "The input vcf file")
argumentDef+=("-database" database "The species database to use. The current options are:\n`ls $snpeffDir/data/ | sed 's:/::'`")
argumentDef+=("-cancerHeader" cancerHead "Indicates that you would like to pass the cancer flag, and that you have added pedigree information to the vcf header" 0 "no")
argumentDef+=("-cancerTextFile" cancerText "Indicates that you would like to pass the cancer flag, and indicates the text file which lists the germline / somatic samples" 1 "")
argumentDef+=("-odir" outDir "The place to put the output files" 1 "./")

argumentDef+=("-name" geName "Name of the submitted grid engine job" 1 "snpEff")
argumentDef+=("-email" emailA "Email address to be emailed when the job completes." 1 "")

argumentsArr=("$@")
handleArguments "$description" argumentDef[@] argumentsArr[@]

convertPathToAbsolute $vcfFile vcfFile
convertPathToAbsolute $outDir outDir

if [ "$emailA" != "" ]; then
    emailA="-m ea -M $emailA"
fi

outbase=`basename $vcfFile | sed 's:\.[^\.]*$::'`
cancerArg=""
if [ "$cancerHead" == "yes" ]; then
    if [ "$cancerText" == "yes" ]; then
	echo "Please pass only one of \"cancerHeader\" and \"cancerTextFile\""
	exit 1
    fi
    cancerArg="-cancer"
elif [ "$cancerText" != "" ]; then
    convertPathToAbsolute $cancerText cancerText
    cancerArg="-cancer -cancerSamples $cancerText"
fi
cd $outDir
qsub -cwd -b y -V -l h_vmem=5G -N $geName -e $outDir/$outbase.err -o $outDir/$outbase.snpeff.vcf $emailA java -Xmx4g -jar $snpeffDir/snpEff.jar eff -c $snpeffDir/snpEff.config -v $cancerArg $database $vcfFile -formatEff 
