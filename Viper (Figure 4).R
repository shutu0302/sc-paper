#subset
cHCs3_IHCs_subset <- subset(combined6, ident=c("cHCs3","IHCs"))
#VCT
cHCs3_IHCs_subset <- SCTransform(cHCs3_IHCs_subset, verbose = FALSE)
#GET assay data
c3_IHC_matrix <- GetAssayData(object = cHCs3_IHCs_subset, slot = "scale.data")

saveRDS(two_matrix, file = "two_matrix.Rds")
saveRDS(centroids, file = "centroids.Rds")



###viper analysis

signature <- rowTtest(all_ed, "cell.type", c("cHCs1","cHCs2","cHCs3"), c("IHCs","OHCs"))
signature <- (qnorm(signature$p.value/2, lower.tail = FALSE),sign(signature$statistic))[, 1]
nullmodel <- ttestNull(all_ed, "cell.type", c("cHCs1","cHCs2","cHCs3"), c("IHCs","OHCs"), per = 1000,repos = TRUE, verbose = FALSE)
mrs <- msviper(signature, regulon_for_allcell, nullmodel, verbose = FALSE)
plot(mrs)
plot(mrs,mrs=20)


mrs <- ledge(mrs)
signature <- bootstrapTtest(all_ed, "cell.type", c("cHCs1","cHCs2","cHCs3"), c("IHCs","OHCs"), verbose = FALSE)
mrs <- msviper(signature, regulon_for_allcell, nullmodel, verbose = FALSE)
mrs <- bootstrapmsviper(mrs, "mode")
plot(mrs, cex = .7)
plot(mrs,mrs=20)

#synergistic effect
mrs <- msviperCombinatorial(mrs, regulators = 63, verbose = FALSE)
mrs <- msviperSynergy(mrs, verbose = FALSE)
summary(mrs)
