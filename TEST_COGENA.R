library(cogena)

data(Psoriasis)

# KEGG Pathway gene set
annoGMT <- "c2.cp.kegg.v7.01.symbols.gmt.xz"
# GO biological process gene set
# annoGMT <- "c5.bp.v7.0.symbols.gmt.xz"
annofile <- system.file("extdata", annoGMT, package="cogena")
# the number of clusters. It can be a vector.
# nClust <- 2:20
nClust <- 10
# Making factor "Psoriasis" behind factor "ct" means Psoriasis Vs Control
# is up-regualted
sampleLabel <- factor(sampleLabel, levels=c("ct", "Psoriasis"))
# the number of cores.
# ncore <- 8
ncore <- 6
# the clustering methods
# clMethods <- c("hierarchical","kmeans","diana","fanny","som","model",
# "sota","pam","clara","agnes") # All the methods can be used together.
clMethods <- c("hierarchical","pam")
# the distance metric
metric <- "correlation"
# the agglomeration method used for hierarchical clustering
# (hierarchical and agnes)
method <- "complete"

# Co-expression Analysis
genecl_result <- coExp(DEexprs, nClust=nClust, clMethods=clMethods,
                       metric=metric, method=method, ncore=ncore)

# Enrichment (Pathway) analysis for the co-expressed genes
clen_res <- clEnrich(genecl_result, annofile=annofile, sampleLabel=sampleLabel)


# Here we consider the "pam" method and 10 clusters.
# Always make the number as character, please!
enrichment.table <- enrichment(clen_res, "pam", "10")


heatmapCluster(clen_res, "pam", "10", maintitle="Psoriasis")


heatmapPEI(clen_res, "pam", "10", printGS=FALSE, maintitle="Pathway analysis for Psoriasis")


###Drug repositioning

# A comprehensive way
# cmapDn100_cogena_result <- clEnrich(genecl_result,
# annofile=system.file("extdata", "CmapDn100.gmt.xz", package="cogena"),
# sampleLabel=sampleLabel)
# A quick way
# Based on the pathway analysis results
cmapDn100_cogena_result <- clEnrich_one(genecl_result, "pam", "10",
                                        annofile=system.file("extdata", "CmapDn100.gmt.xz", package="cogena"),
                                        sampleLabel=sampleLabel)

heatmapPEI(cmapDn100_cogena_result, "pam", "10", printGS=FALSE,
           orderMethod = "7", maintitle="Drug repositioning for Psoriasis")

###MultiInstance
heatmapCmap(cmapDn100_cogena_result, "pam", "10", printGS=FALSE,
            orderMethod = "7", maintitle="Drug repositioning for Psoriasis")
