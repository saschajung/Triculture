###############################################################################
# This script generates the differential expression analysis results of
# OPCs, astrocytes, radial glia and neurons between the Linnarsson data and the triculture.
# 
# Results are shown in Supplementary Table 4.4 of the manuscript
###############################################################################

library(Seurat)
library(WebGestaltR)

load("./IntegrationWithOtherStudies/Linnarsson/DifferentialExpressionResults.RData")
combined.sct <- readRDS("./IntegrationWithOtherStudies/Linnarsson/Combined_SeuratObj.rds")

inp_microglia_list <- list("This study" = microglia.diff$gene[microglia.diff$p_val_adj <= 0.01 & microglia.diff$avg_logFC > 0.5],
                           "Brain" = microglia.diff$gene[microglia.diff$p_val_adj <= 0.01 & microglia.diff$avg_logFC < -0.5])

micBrain <- WebGestaltR(enrichMethod = "ORA",
                        organism = "hsapiens",
                        enrichDatabase = "geneontology_Biological_Process",
                        interestGene = inp_microglia_list$Brain,
                        interestGeneType = "ensembl_gene_id",
                        referenceGene = rownames(combined.sct@assays$RNA),
                        referenceGeneType = "ensembl_gene_id",
                        sigMethod = "fdr",
                        fdrMethod = "BH",
                        fdrThr = 0.01,
                        isOutput = F
)

micThis <- WebGestaltR(enrichMethod = "ORA",
                       organism = "hsapiens",
                       enrichDatabase = "geneontology_Biological_Process",
                       interestGene = inp_microglia_list$`This study`,
                       interestGeneType = "ensembl_gene_id",
                       referenceGene = rownames(combined.sct@assays$RNA),
                       referenceGeneType = "ensembl_gene_id",
                       sigMethod = "fdr",
                       fdrMethod = "BH",
                       fdrThr = 0.01,
                       isOutput = F
)

#Excitatory Neuron
inp_excitneuron_list <- list("This study" = excitNeuron.diff$gene[excitNeuron.diff$p_val_adj <= 0.01 & excitNeuron.diff$avg_logFC > 0.5],
                             "Brain" = excitNeuron.diff$gene[excitNeuron.diff$p_val_adj <= 0.01 & excitNeuron.diff$avg_logFC < -0.5])

excitNeuBrain <- WebGestaltR(enrichMethod = "ORA",
                             organism = "hsapiens",
                             enrichDatabase = "geneontology_Biological_Process",
                             interestGene = inp_excitneuron_list$Brain,
                             interestGeneType = "ensembl_gene_id",
                             referenceGene = rownames(combined.sct@assays$RNA),
                             referenceGeneType = "ensembl_gene_id",
                             sigMethod = "fdr",
                             fdrMethod = "BH",
                             fdrThr = 0.01,
                             isOutput = F
)

excitNeuThis <- WebGestaltR(enrichMethod = "ORA",
                            organism = "hsapiens",
                            enrichDatabase = "geneontology_Biological_Process",
                            interestGene = inp_excitneuron_list$`This study`,
                            interestGeneType = "ensembl_gene_id",
                            referenceGene = rownames(combined.sct@assays$RNA),
                            referenceGeneType = "ensembl_gene_id",
                            sigMethod = "fdr",
                            fdrMethod = "BH",
                            fdrThr = 0.01,
                            isOutput = F
)

#Inhibitory Neuron
inp_inhibneuron_list <- list("This study" = inhibNeuron.diff$gene[inhibNeuron.diff$p_val_adj <= 0.01 & inhibNeuron.diff$avg_logFC > 0.5],
                             "Brain" = inhibNeuron.diff$gene[inhibNeuron.diff$p_val_adj <= 0.01 & inhibNeuron.diff$avg_logFC < -0.5])

inhibNeuBrain <- WebGestaltR(enrichMethod = "ORA",
                             organism = "hsapiens",
                             enrichDatabase = "geneontology_Biological_Process",
                             interestGene = inp_inhibneuron_list$Brain,
                             interestGeneType = "ensembl_gene_id",
                             referenceGene = rownames(combined.sct@assays$RNA),
                             referenceGeneType = "ensembl_gene_id",
                             sigMethod = "fdr",
                             fdrMethod = "BH",
                             fdrThr = 0.01,
                             isOutput = F
)

inhibNeuThis <- WebGestaltR(enrichMethod = "ORA",
                            organism = "hsapiens",
                            enrichDatabase = "geneontology_Biological_Process",
                            interestGene = inp_inhibneuron_list$`This study`,
                            interestGeneType = "ensembl_gene_id",
                            referenceGene = rownames(combined.sct@assays$RNA),
                            referenceGeneType = "ensembl_gene_id",
                            sigMethod = "fdr",
                            fdrMethod = "BH",
                            fdrThr = 0.01,
                            isOutput = F
)

#Astrocyte
inp_astro_list <- list("This study" = astrocyte.diff$gene[astrocyte.diff$p_val_adj <= 0.01 & astrocyte.diff$avg_logFC > 0.5],
                       "Brain" = astrocyte.diff$gene[astrocyte.diff$p_val_adj <= 0.01 & astrocyte.diff$avg_logFC < -0.5])

astroBrain <- WebGestaltR(enrichMethod = "ORA",
                          organism = "hsapiens",
                          enrichDatabase = "geneontology_Biological_Process",
                          interestGene = inp_astro_list$Brain,
                          interestGeneType = "ensembl_gene_id",
                          referenceGene = rownames(combined.sct@assays$RNA),
                          referenceGeneType = "ensembl_gene_id",
                          sigMethod = "fdr",
                          fdrMethod = "BH",
                          fdrThr = 0.01,
                          isOutput = F
)

astroThis <- WebGestaltR(enrichMethod = "ORA",
                         organism = "hsapiens",
                         enrichDatabase = "geneontology_Biological_Process",
                         interestGene = inp_astro_list$`This study`,
                         interestGeneType = "ensembl_gene_id",
                         referenceGene = rownames(combined.sct@assays$RNA),
                         referenceGeneType = "ensembl_gene_id",
                         sigMethod = "fdr",
                         fdrMethod = "BH",
                         fdrThr = 0.01,
                         isOutput = F
)

#OPC
inp_opc_list <- list("This study" = opc.diff$gene[opc.diff$p_val_adj <= 0.01 & opc.diff$avg_logFC > 0.5],
                     "Brain" = opc.diff$gene[opc.diff$p_val_adj <= 0.01 & opc.diff$avg_logFC < -0.5])

opcBrain <- WebGestaltR(enrichMethod = "ORA",
                        organism = "hsapiens",
                        enrichDatabase = "geneontology_Biological_Process",
                        interestGene = inp_opc_list$Brain,
                        interestGeneType = "ensembl_gene_id",
                        referenceGene = rownames(combined.sct@assays$RNA),
                        referenceGeneType = "ensembl_gene_id",
                        sigMethod = "fdr",
                        fdrMethod = "BH",
                        fdrThr = 0.01,
                        isOutput = F
)

opcThis <- WebGestaltR(enrichMethod = "ORA",
                       organism = "hsapiens",
                       enrichDatabase = "geneontology_Biological_Process",
                       interestGene = inp_opc_list$`This study`,
                       interestGeneType = "ensembl_gene_id",
                       referenceGene = rownames(combined.sct@assays$RNA),
                       referenceGeneType = "ensembl_gene_id",
                       sigMethod = "fdr",
                       fdrMethod = "BH",
                       fdrThr = 0.01,
                       isOutput = F
)

save(astroBrain,astroThis,micBrain,micThis,excitNeuBrain,excitNeuThis,inhibNeuBrain,inhibNeuThis,opcThis,opcBrain,file = "./IntegrationWithOtherStudies/Linnarsson/GOResults.RData")

write.table(astroBrain,file = "./IntegrationWithOtherStudies/Linnarsson/GO_AstrocyteBrain.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(astroThis,file = "./IntegrationWithOtherStudies/Linnarsson/GO_AstrocyteTriculture.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(micBrain,file = "./IntegrationWithOtherStudies/Linnarsson/GO_MicrogliaBrain.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(micThis,file = "./IntegrationWithOtherStudies/Linnarsson/GO_MicrogliaTriculture.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(excitNeuBrain,file = "./IntegrationWithOtherStudies/Linnarsson/GO_ExcitNeuBrain.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(excitNeuThis,file = "./IntegrationWithOtherStudies/Linnarsson/GO_ExcitNeuTriculture.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(inhibNeuBrain,file = "./IntegrationWithOtherStudies/Linnarsson/GO_InhibNeuBrain.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(inhibNeuThis,file = "./IntegrationWithOtherStudies/Linnarsson/GO_InhibNeuTriculture.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(opcBrain,file = "./IntegrationWithOtherStudies/Linnarsson/GO_OPCBrain.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(opcThis,file = "./IntegrationWithOtherStudies/Linnarsson/GO_OPCTriculture.txt",sep = "\t",row.names = F, col.names = T, quote = F)
