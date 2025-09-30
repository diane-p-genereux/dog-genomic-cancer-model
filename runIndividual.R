#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

#plain text name for files
name <- args[1]
#variant file, format Sample / Ref / Alt / Trinuc context (with ref allele)
variants <- args[2]
#number of samples
#Signature model to use; default is "multinomial"
model <- "multinomial"
#Where to put the files
out_dir <- args[3]
numIter <- as.integer(args[4])

library(sigfit)

data(cosmic_signatures_v3)

# Load variant data and make mutational matrix
currFile <- read.table(variants, header=T, sep="\t", colClasses=c("character", "character", "character"))
varMatrix = build_catalogues(currFile[currFile[,2] != currFile[,3],])
#trinuc R files that I made; details in /seq/vgb/jason/cancer_mutations/cancer_machine_learning/cancer_unsupervised/signatures/
# Read in trinucleotide contexts for genome in question
if(grepl("canFam3", variants)) {
	oppsMatrix <- get(load("/seq/vgb/jason/cancer_mutations/cancer_machine_learning/cancer_unsupervised/original_attempt_with_FAMD/signatures/trinuc_freqs_canFam3.1.RData"))
} else if(grepl("hg38", variants)) {
	oppsMatrix <- get(load("/seq/vgb/jason/cancer_mutations/cancer_machine_learning/cancer_unsupervised/original_attempt_with_FAMD/signatures/trinuc_freqs_hg38.RData"))
} else if(grepl("hg19", variants)) {
	oppsMatrix <- get(load("/seq/vgb/jason/cancer_mutations/cancer_machine_learning/cancer_unsupervised/original_attempt_with_FAMD/signatures/trinuc_freqs_hg19.RData"))
} else {
	stop(paste("unidentified genome", variants))
}

# Normalize the human signatures to be genome-agnostic
cosmic.normalised <- convert_signatures(cosmic_signatures_v3, opportunities_from="human-genome")

# Fit signatures
mcmc_samples_fit <- fit_signatures(counts = varMatrix, opportunities=oppsMatrix, signatures = cosmic.normalised, model = model, iter = numIter, warmup = numIter/2, chains = 1, seed = 1756)

# Get exposures
exposures <- retrieve_pars(mcmc_samples_fit, par = "exposures", hpd_prob = 0.95)


# Get lower end of Bayesian HPD interval of exposures
lower <- exposures$lower_95
write.table(lower, file=paste0(out_dir, "lower/", name,"_lower.txt", collapse=NULL), sep="\t", col.names=NA, quote=FALSE)

selected <- apply(lower, 2, function (X) any(X >= 0.01))
write.table(selected, file=paste0(out_dir, "selected/", name, "_selectedSigs.txt", collapse=NULL), sep="\t", col.names=FALSE, quote=FALSE)
if(length(selected[selected]) == 1) {
	fileConn<-file(paste0(out_dir, "exposures/", name, "_exposures.txt", collapse=NULL))
	writeLines(c(paste0("\t", names(selected[selected])), "1.0"), fileConn)
	close(fileConn)
} else if(length(selected[selected]) > 1) {
	if(any(rownames(selected) != rownames(cosmic.normalised))) {
		stop("rownames don't match")
	}
	cosmic.normalised <- cosmic.normalised[unlist(selected), ]

	mcmc_samples_fit <- fit_signatures(counts = varMatrix, opportunities=oppsMatrix, signatures = cosmic.normalised, model = model, iter = numIter, warmup = numIter/2, chains = 1, seed = 1756)

	# Get exposures
	exposures <- retrieve_pars(mcmc_samples_fit, par = "exposures", hpd_prob = 0.95)

	write.table(exposures$mean, file=paste0(out_dir, "exposures/", name,"_exposures.txt", collapse=NULL), sep="\t", col.names=NA, quote=FALSE)
	plot_all(mcmc_samples = mcmc_samples_fit, out_path = paste0(out_dir, "pdfs/"), prefix = name)
} else {
	fileConn<-file(paste0(out_dir, "exposures/", name, "_exposures.txt", collapse=NULL))
	writeLines(c("No significant exposures found"), fileConn)
	close(fileConn)
}
