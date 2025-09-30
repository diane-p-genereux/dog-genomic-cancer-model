# dog-genomic-cancer-model

Genomics-guided therapies have transformed clinical outcomes for some cancer types, but development of new treatments remains slow. Pet dogs are an underutilized and potentially powerful model for therapeutic innovation. Here, we systematically evaluate genomic similarity between dog and human tumors using the largest comparative cancer dataset assembled to date: 15,315 orthologous genes in 429 dog tumors and 14,966 human tumors across 39 different cancer types. We find that cancers in dogs and humans are genomically almost indistinguishable, with shared mutational signatures and recurrent mutations in the same genes. Some cancer driver genes are more frequently mutated in dogs than in humans, creating an opportunity to recruit large cohorts for clinical trials. Tumors from dogs and humans are, on average, as genomically similar as tumors from different human cancer types, and tumors do not separate by species in unsupervised clustering. Supervised machine learning identified shared genomic features between dog and human tumors. Even so, consistent with the genomic heterogeneity of both dog and human cancers, there was no dog cancer type for which all tumors were classified to a single human type. Our findings establish dogs as a compelling system for the development of targeted cancer therapies and illustrate the potential for using machine learning to discover models for human cancers in dogs and other species. 



Code for genomic comparisons of dog to human cancers, as implemented in Genereux DP, Megquier K, et al., and Karlsson EK.



Files: 

1. Parse variants from Variant Call Format, Mutation Annotation Format (MAF), or custom formats using a custom python script leveraging the pysam14 wrapper of the Samtools15,16 and htslib17 tools.
   parseVCF_110723.py
3. Run the SNPEff variant-effect predictor on the dog cancer samples.: 
    runSNPEFF*.sh
   snpeff_ge_NEW_TERMS_FOR_ML.sh
4. Calculate dN/dS scores to detect genes under selection:
     dndscvExample.R
6. Screen samples for mutational signatures
   cancer_machine_learning/mutational_signatures/runIndividual.R
7. Apply machine-learning models for feature selection model training: 
  chooseInformativeFeaturesAndTrainModelsUpdate.py
8. Implement machine learing classification and calculate Shapley values:
   runModelUpdate2.py
