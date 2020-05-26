library(MutationalPatterns)
library(GenomicRanges)
library(VariantAnnotation)
library(SomaticSignatures)
#library(deconstructSigs)

###################################################################
# Load lacZ-corrected COSMIC signatures (V3) WITHOUT control
lacZ_cosmic_signatures_v3_NO_control <- as.data.frame(t(read.table("data/lacZ_corrected_COSMIC_signatures_v3_without_control.txt",
                                                                sep="\t",
                                                                header=F,
                                                                stringsAsFactors = F)))
newSigs <- read.table("data/ColumnNamesforSignatures3withoutcontrol.txt", header=T)
rownames(lacZ_cosmic_signatures_v3_NO_control) <- colnames(newSigs)
colnames(lacZ_cosmic_signatures_v3_NO_control) <- colnames(signatures.cosmic)

###################################################################
# Load lacZ-corrected COSMIC signatures (V3) with control
lacZ_cosmic_signatures_v3_control <- as.data.frame(t(read.table("data/lacZ_corrected_COSMIC_signatures_v3_with_control.txt",
                                                     sep="\t",
                                                     header=T,
                                                     stringsAsFactors = F)))
colnames(lacZ_cosmic_signatures_v3_control) <- colnames(signatures.cosmic)

###################################################################
# Load lacZ-corrected COSMIC signatures (V2) with control
lacZ_cosmic_signatures_v2_control <- as.data.frame(t(read.table("data/lacZ_corrected_COSMIC_signatures_v2_with_control.txt",
                                                                     sep="\t",
                                                                     header=T,
                                                                     stringsAsFactors = F)))
colnames(lacZ_cosmic_signatures_v2_control) <- colnames(signatures.cosmic)

###################################################################
# Load mutation observations
mutationData <- read.table("./Meta Data.txt",
                           header=T,
                           sep="\t")
mutationData$chr <- "lacZ"
colnames(mutationData)[2] <- "start"
mutationData$end <- mutationData$start

mutationDataVRange <- VRanges(seqnames = mutationData$chr,
                                ranges = IRanges(mutationData$start,
                                                 end = mutationData$end),
                                ref = mutationData$Ref,
                                alt = mutationData$Alt,
                                sampleNames = mutationData$Group)
###################################################################
# Generate motif matrix for mutation data
lacZ.fa = FaFile(paste("data/lacZ.fa", sep=""))
mutationRanges <- mutationContext(mutationDataVRange, lacZ.fa)
mut_matrix = as.matrix(motifMatrix(mutationRanges))

###################################################################
# Do fitting of signatures
fit_res_v2 <- fit_to_signatures(mut_matrix, t(lacZ_cosmic_signatures_v2_control))
fit_res_v3 <- fit_to_signatures(mut_matrix, t(lacZ_cosmic_signatures_v3_control))
fit_res_v3_no_control <- fit_to_signatures(mut_matrix, t(lacZ_cosmic_signatures_v3_NO_control))

###################################################################
# Plot results
cutoff=0.0

# V2 Signatures, with control signature
selectv2 <- which(rowSums(fit_res_v2$contribution) > cutoff)
plot_contribution(fit_res_v2$contribution[selectv2,],
                  lacZ_cosmic_signatures_v2_control[,selectv2],
                  coord_flip = FALSE,
                  mode = "relative")
plot_contribution_heatmap(fit_res_v2$contribution,
                          cluster_samples = TRUE,
                          method = "complete")


# V3 Signatures, with control signature
selectv3 <- which(rowSums(fit_res_v3$contribution) > cutoff)
plot_contribution(fit_res_v3$contribution[selectv3,],
                  lacZ_cosmic_signatures_v3_control[,selectv3],
                  coord_flip = FALSE,
                  mode = "relative")
plot_contribution_heatmap(fit_res_v3$contribution,
                          cluster_samples = TRUE,
                          method = "complete")


# V3 Signatures, without control signature
selectv3_no_control <- which(rowSums(fit_res_v3_no_control$contribution) > cutoff)
plot_contribution(fit_res_v3_no_control$contribution[selectv3_no_control,],
                  lacZ_cosmic_signatures_v3_NO_control[,selectv3_no_control],
                  coord_flip = FALSE,
                  mode = "relative")
plot_contribution_heatmap(fit_res_v3_no_control$contribution,
                          cluster_samples = TRUE,
                          method = "complete")

for (i in 1:ncol(mut_matrix)) {
print(plot_compare_profiles(mut_matrix[,i], fit_res_v2$reconstructed[,i],
                      profile_names = c("Original", "Reconstructed"),
                      condensed = TRUE))
plot_compare_profiles(mut_matrix[,i], fit_res_v3$reconstructed[,i],
                      profile_names = c("Original", "Reconstructed"),
                      condensed = TRUE)
plot_compare_profiles(mut_matrix[,i], fit_res_v3_no_control$reconstructed[,i],
                      profile_names = c("Original", "Reconstructed"),
                      condensed = TRUE)
}
