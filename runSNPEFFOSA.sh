#!/bin/bash                                                                                                                                                                                                                  
source /broad/software/scripts/useuse
use Java-1.8
use UGER

mkdir SNPeff_output

for vcfFile in /seq/vgb/sakthi/b-OSA_WES_TN_UCSC_uploads/m2pass_VCF/IndivVCF/*.vcf
do    
    database=Canis_lupus_familiaris
    odir=/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/OSA/
    id=`basename $vcfFile | sed 's:.vcf::'`
    if [ ! -f ./OSA/$id.snpeff.vcf ]; then
	echo $id
	/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/snpeff_ge_NEW_TERMS_FOR_ML.sh -vcf $vcfFile -database $database -cancerTextFile -odir $odir
    fi
done

