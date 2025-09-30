#!/bin/bash                                                                                                                                                                                                                  
source /broad/software/scripts/useuse
use Java-1.8
use UGER

#mkdir SNPeff_output

#cd ./LSA/temp/
#/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/LSATest/*.tgz
#/seq/vgb/jason/cancer_mutations/ginger_lymphoma_somatic_mutations/sample_runs/*/*.snps.tgz
#cd ./LSATest/temp/

#for file in /seq/vgb/jason/cancer_mutations/ginger_lymphoma_somatic_mutations/sample_runs/*/*.snps.tgz
#do
#    id=$(basename $file .snps.tgz)
#    if [ ! -f ./$id.snp.filtered.vcf ]; then
	#echo $id.snp.filtered.vcf
#	tar -C ./ -xvf $file
#	cp ./$id.snps/*.filtered.vcf ./
#	rm -r ./$id.snps
#    fi
#done




for vcfFile in /seq/vgb/cancer_r01/maja_cmt/recmtdata/FINAL_VCF_PON_IDOG_ERIK_FILT/*vcf.gz
do
    database=Canis_lupus_familiaris
    odir=/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/CMT/
    /seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/snpeff_ge_NEW_TERMS_FOR_ML.sh -vcf $vcfFile -database $database -cancerTextFile -odir $odir
done

