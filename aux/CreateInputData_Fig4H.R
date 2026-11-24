
if(!file.exists("./IntegrationWithOtherStudies/scOrganoids/GOEnrichment.RData")){
  source("./Analyses/GOEnrichment_w_scOrganoids.R")
}

load("./IntegrationWithOtherStudies/scOrganoids/GOEnrichment.RData",verbose = T)

syngo <- read.table("./InputData/syngo_ontologies_BP_GO.txt",header = F)
neuronBrain_SynGO <- neuronBrain[which(neuronBrain$geneSet %in% syngo$V1),]
neuronThis_SynGO <- neuronThis[which(neuronThis$geneSet %in% syngo$V1),]

saveRDS(neuronThis_SynGO,"./InputData/InputData_Fig4H.rds")
