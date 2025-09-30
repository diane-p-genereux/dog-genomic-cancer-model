# dog-genomic-cancer-model
Code for genomic comparisons of dog to human cancers, as implemented in Genereux DP, Megquier K, et al., and Karlsson EK.



Files: 

1. Code to parse variants from Variant Call Format, Mutation Annotation Format (MAF), or custom formats using a custom python script leveraging the pysam14 wrapper of the Samtools15,16 and htslib17 tools.
   parseVCF_110723.py
3. Code torun the SNPEff variant-effect predictor on the dog cancer samples.: 
    runSNPEFF*.sh
   snpeff_ge_NEW_TERMS_FOR_ML.sh
4. Calculate dN/dS scores to detect genes under selection:
     dndscvExample.R
6. Screens samples for mutational signatures
   cancer_machine_learning/mutational_signatures/runIndividual.R
7. Apply machine-learning models for feature selection model training: 
  chooseInformativeFeaturesAndTrainModelsUpdate.py
8. Implementing machine learing classification and calculate Shapley values:
   runModelUpdate2.py
