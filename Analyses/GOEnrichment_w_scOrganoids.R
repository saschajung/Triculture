###############################################################################
# This script generates the GO enrichment results of
# astrocytes, radial glia and neurons both in organoids and the triculture.
# 
# Results are shown in Supplementary Table 4.6 of the manuscript
###############################################################################

library(WebGestaltR)

if(!file.exists("./IntegrationWithOtherStudies/scOrganoids/DiffExp.RData")){
  source("./Analyses/DiffExp_w_scOrganoids.R")
}

load("./IntegrationWithOtherStudies/scOrganoids/DiffExp.RData",verbose = T)

# Neuron
inp_neuron_list <- list("This study" = neuron.diff$gene[neuron.diff$p_val_adj <= 0.01 & neuron.diff$avg_logFC > 0.5],
                        "Brain" = neuron.diff$gene[neuron.diff$p_val_adj <= 0.01 & neuron.diff$avg_logFC < -0.5])

neuronBrain <- WebGestaltR(enrichMethod = "ORA",
                           organism = "hsapiens",
                           enrichDatabase = "geneontology_Biological_Process",
                           interestGene = inp_neuron_list$Brain,
                           interestGeneType = "genesymbol",
                           referenceGene = neuron.diff$gene,
                           referenceGeneType = "genesymbol",
                           sigMethod = "fdr",
                           fdrMethod = "BH",
                           fdrThr = 0.01,
                           isOutput = F
)

neuronThis <- WebGestaltR(enrichMethod = "ORA",
                          organism = "hsapiens",
                          enrichDatabase = "geneontology_Biological_Process",
                          interestGene = inp_neuron_list$`This study`,
                          interestGeneType = "genesymbol",
                          referenceGene = neuron.diff$gene,
                          referenceGeneType = "genesymbol",
                          sigMethod = "fdr",
                          fdrMethod = "BH",
                          fdrThr = 0.01,
                          isOutput = F
)

# Astrocyte
inp_astro_list <- list("This study" = astrocyte.diff$gene[astrocyte.diff$p_val_adj <= 0.01 & astrocyte.diff$avg_logFC > 0.5],
                       "Brain" = astrocyte.diff$gene[astrocyte.diff$p_val_adj <= 0.01 & astrocyte.diff$avg_logFC < -0.5])

astroBrain <- WebGestaltR(enrichMethod = "ORA",
                          organism = "hsapiens",
                          enrichDatabase = "geneontology_Biological_Process",
                          interestGene = inp_astro_list$Brain,
                          interestGeneType = "genesymbol",
                          referenceGene = astrocyte.diff$gene,
                          referenceGeneType = "genesymbol",
                          sigMethod = "fdr",
                          fdrMethod = "BH",
                          fdrThr = 0.01,
                          isOutput = F
)

astroThis <- WebGestaltR(enrichMethod = "ORA",
                         organism = "hsapiens",
                         enrichDatabase = "geneontology_Biological_Process",
                         interestGene = inp_astro_list$`This study`,
                         interestGeneType = "genesymbol",
                         referenceGene = astrocyte.diff$gene,
                         referenceGeneType = "genesymbol",
                         sigMethod = "fdr",
                         fdrMethod = "BH",
                         fdrThr = 0.01,
                         isOutput = F
)

# RG
inp_rg_list <- list("This study" = radialglia.diff$gene[radialglia.diff$p_val_adj <= 0.01 & radialglia.diff$avg_logFC > 0.5],
                    "Brain" = radialglia.diff$gene[radialglia.diff$p_val_adj <= 0.01 & radialglia.diff$avg_logFC < -0.5])

rgBrain <- WebGestaltR(enrichMethod = "ORA",
                       organism = "hsapiens",
                       enrichDatabase = "geneontology_Biological_Process",
                       interestGene = inp_rg_list$Brain,
                       interestGeneType = "genesymbol",
                       referenceGene = radialglia.diff$gene,
                       referenceGeneType = "genesymbol",
                       sigMethod = "fdr",
                       fdrMethod = "BH",
                       fdrThr = 0.01,
                       isOutput = F
)

rgThis <- WebGestaltR(enrichMethod = "ORA",
                      organism = "hsapiens",
                      enrichDatabase = "geneontology_Biological_Process",
                      interestGene = inp_rg_list$`This study`,
                      interestGeneType = "genesymbol",
                      referenceGene = radialglia.diff$gene,
                      referenceGeneType = "genesymbol",
                      sigMethod = "fdr",
                      fdrMethod = "BH",
                      fdrThr = 0.01,
                      isOutput = F
)

save(astroBrain,astroThis,rgBrain,rgThis,neuronBrain,neuronThis,file = "./IntegrationWithOtherStudies/scOrganoids/GOEnrichment.RData")

write.table(astroBrain,file = "./IntegrationWithOtherStudies/scOrganoids/GO_AstrocyteOrganoid.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(astroThis,file = "./IntegrationWithOtherStudies/scOrganoids/GO_AstrocyteTriculture.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(rgBrain,file = "./IntegrationWithOtherStudies/scOrganoids/GO_RadialGliaOrganoid.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(rgThis,file = "./IntegrationWithOtherStudies/scOrganoids/GO_RadialGliaTriculture.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(neuronBrain,file = "./IntegrationWithOtherStudies/scOrganoids/GO_NeuronOrganoid.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(neuronThis,file = "./IntegrationWithOtherStudies/scOrganoids/GO_NeuronTriculture.txt",sep = "\t",row.names = F, col.names = T, quote = F)

