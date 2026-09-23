#ORGANOIDS

library(Libra)
load("./OtherStudies/scOrganoids/Y1Integration/Combined_SeuratObj.RData",verbose = T)
combined.sct$Origin <- "This study"
combined.sct$Origin[combined.sct$orig.ident != "Y1"] <- "Organoid"

combined.sct$CTOrigin <- paste0(combined.sct$CTs," (",combined.sct$Origin,")")

combined.sct$Clusts <- paste0(combined.sct$seurat_clusters," (",combined.sct$Origin,")")

sub <- subset(combined.sct, subset = CTOrigin %in% c("Astrocyte (This study)","Astroglia (Organoid)"))
DefaultAssay(sub) <- "RNA"
sub$CTs <- "Astrocyte"
sub$replicate <- sub$orig.ident
set.seed(1234)
sub$replicate[which(sub$orig.ident == "11a.6mon")] <- sample(paste0("11a.6mon",1:2),length(which(sub$orig.ident == "11a.6mon")),replace = T)
sub$replicate[which(sub$orig.ident == "PGP1.6mon")] <- sample(paste0("PGP1.6mon",1:2),length(which(sub$orig.ident == "PGP1.6mon")),replace = T)
sub$replicate[which(sub$orig.ident == "Y1")] <- sample(paste0("Y1",1:4),length(which(sub$orig.ident == "Y1")),replace = T)
for(assay_name in Assays(sub)) {
  sub[[assay_name]] <- as(sub[[assay_name]], "Assay5")
}
astrocyte.diff <- Libra::run_de(input = sub,
                                cell_type_col = "CTs",
                                replicate_col = "replicate",
                                label_col = "Origin",
                                de_family = "pseudobulk",
                                de_method = "edgeR",
                                de_type = "QLF") #Negative FC means higher in organoid

sub <- subset(combined.sct, subset = CTOrigin %in% c("Neuron (This study)","CPNs (Organoid)"))
DefaultAssay(sub) <- "RNA"
sub$CTs <- "Neuron"
sub$replicate <- sub$orig.ident
set.seed(1234)
sub$replicate[which(sub$orig.ident == "11a.6mon")] <- sample(paste0("11a.6mon",1:2),length(which(sub$orig.ident == "11a.6mon")),replace = T)
sub$replicate[which(sub$orig.ident == "PGP1.6mon")] <- sample(paste0("PGP1.6mon",1:2),length(which(sub$orig.ident == "PGP1.6mon")),replace = T)
sub$replicate[which(sub$orig.ident == "Y1")] <- sample(paste0("Y1",1:4),length(which(sub$orig.ident == "Y1")),replace = T)
for(assay_name in Assays(sub)) {
  sub[[assay_name]] <- as(sub[[assay_name]], "Assay5")
}
neuron.diff <- Libra::run_de(input = sub,
                             cell_type_col = "CTs",
                             replicate_col = "replicate",
                             label_col = "Origin",
                             de_family = "pseudobulk",
                             de_method = "edgeR",
                             de_type = "QLF") #Negative FC means higher in organoid

sub <- subset(combined.sct, subset = CTOrigin %in% c("Radial Glia (This study)","RG (Organoid)"))
DefaultAssay(sub) <- "RNA"
sub$CTs <- "Radial Glia"
sub$replicate <- sub$orig.ident
set.seed(1234)
sub$replicate[which(sub$orig.ident == "11a.6mon")] <- sample(paste0("11a.6mon",1:2),length(which(sub$orig.ident == "11a.6mon")),replace = T)
sub$replicate[which(sub$orig.ident == "PGP1.6mon")] <- sample(paste0("PGP1.6mon",1:2),length(which(sub$orig.ident == "PGP1.6mon")),replace = T)
sub$replicate[which(sub$orig.ident == "Y1")] <- sample(paste0("Y1",1:4),length(which(sub$orig.ident == "Y1")),replace = T)
for(assay_name in Assays(sub)) {
  sub[[assay_name]] <- as(sub[[assay_name]], "Assay5")
}
radialglia.diff <- Libra::run_de(input = sub,
                                 cell_type_col = "CTs",
                                 replicate_col = "replicate",
                                 label_col = "Origin",
                                 de_family = "pseudobulk",
                                 de_method = "edgeR",
                                 de_type = "QLF") #Negative FC means higher in organoid

sub <- subset(combined.sct, subset = CTOrigin %in% c("Radial Glia (This study)","RG (Organoid)"))


save(astrocyte.diff,radialglia.diff,neuron.diff,file = "./OtherStudies/scOrganoids/Y1Integration/DiffExp_rev.RData",compression_level = 6)
load("./OtherStudies/scOrganoids/Y1Integration/DiffExp_rev.RData", verbose = T)

write.table(astrocyte.diff[which(abs(astrocyte.diff$avg_logFC) > 0.5 & astrocyte.diff$p_val_adj <= 0.01),],"./OtherStudies/scOrganoids/Y1Integration/DiffExp_Astrocyte.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(neuron.diff[which(abs(neuron.diff$avg_logFC) > 0.5 & neuron.diff$p_val_adj <= 0.01),],"./OtherStudies/scOrganoids/Y1Integration/DiffExp_Neuron.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(radialglia.diff[which(abs(radialglia.diff$avg_logFC) > 0.5 & radialglia.diff$p_val_adj <= 0.01),],"./OtherStudies/scOrganoids/Y1Integration/DiffExp_RadialGlia.txt",sep = "\t",row.names = F, col.names = T, quote = F)


library(WebGestaltR)
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

write.table(astroBrain,file = "./OtherStudies/scOrganoids/Y1Integration/GO_AstrocyteOrganoid_rev.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(astroThis,file = "./OtherStudies/scOrganoids/Y1Integration/GO_AstrocyteThis_rev.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(rgBrain,file = "./OtherStudies/scOrganoids/Y1Integration/GO_RadialGliaOrganoid_rev.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(rgThis,file = "./OtherStudies/scOrganoids/Y1Integration/GO_RadialGliaThis_rev.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(neuronBrain,file = "./OtherStudies/scOrganoids/Y1Integration/GO_NeuronOrganoid_rev.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(neuronThis,file = "./OtherStudies/scOrganoids/Y1Integration/GO_NeuronThis_rev.txt",sep = "\t",row.names = F, col.names = T, quote = F)

syngo <- read.table("./syngo_ontologies_BP_GO.txt",header = F)
neuronBrain_SynGO <- neuronBrain[which(neuronBrain$geneSet %in% syngo$V1),]
neuronThis_SynGO <- neuronThis[which(neuronThis$geneSet %in% syngo$V1),]

saveRDS(neuronThis_SynGO,"./PlotsRevision/InputData_Fig4J_rev.rds")

write.table(astroBrain,file = "./OtherStudies/scOrganoids/Y1Integration/AstroBrain.txt",row.names = F,col.names = T,quote =F, sep = "\t")
write.table(astroThis,file = "./OtherStudies/scOrganoids/Y1Integration/AstroThisn.txt",row.names = F,col.names = T,quote =F, sep = "\t")
write.table(neuronBrain,file = "./OtherStudies/scOrganoids/Y1Integration/NeuronBrain.txt",row.names = F,col.names = T,quote =F, sep = "\t")
write.table(neuronThis,file = "./OtherStudies/scOrganoids/Y1Integration/NeuronThis.txt",row.names = F,col.names = T,quote =F, sep = "\t")
write.table(rgBrain,file = "./OtherStudies/scOrganoids/Y1Integration/RGBrain.txt",row.names = F,col.names = T,quote =F, sep = "\t")
write.table(rgThis,file = "./OtherStudies/scOrganoids/Y1Integration/RGThis.txt",row.names = F,col.names = T,quote =F, sep = "\t")

#GSE179590
load("./OtherStudies/GSE179590/sNucIntegration/Combined_SeuratObj.RData",verbose = T)

sub <- subset(combined.sct, subset = CTOrigin %in% c("Microglia (This study)","Microglia (Brain)"))
DefaultAssay(sub) <- "SCT"
sub$replicate <- sub$orig.ident
set.seed(1234)
sub$replicate[which(sub$orig.ident == "snuc")] <- sample(paste0("sNuc",1:2),length(which(sub$orig.ident == "snuc")),replace = T)
for(assay_name in Assays(sub)) {
  sub[[assay_name]] <- as(sub[[assay_name]], "Assay5")
}
microglia.diff <- Libra::run_de(input = sub,
                                meta = sub@meta.data,
                                cell_type_col = "CTs",
                                replicate_col = "replicate",
                                label_col = "Origin",
                                de_family = "pseudobulk",
                                de_method = "edgeR",
                                de_type = "QLF") #Negative FC means higher in this study

sub <- subset(combined.sct, subset = CTOrigin %in% c("Astrocyte (This study)","Astrocyte (Brain)"))
DefaultAssay(sub) <- "SCT"
sub$replicate <- sub$orig.ident
set.seed(1234)
sub$replicate[which(sub$orig.ident == "snuc")] <- sample(paste0("sNuc",1:10),length(which(sub$orig.ident == "snuc")),replace = T)
table(sub$replicate)
for(assay_name in Assays(sub)) {
  sub[[assay_name]] <- as(sub[[assay_name]], "Assay5")
}
astrocyte.diff <- Libra::run_de(input = sub,
                                meta = sub@meta.data,
                                cell_type_col = "CTs",
                                replicate_col = "replicate",
                                label_col = "Origin",
                                de_family = "pseudobulk",
                                de_method = "edgeR",
                                de_type = "QLF") #Negative FC means higher in this study

sub <- subset(combined.sct, subset = CTOrigin %in% c("Neuron (This study)","Excitatory Neuron (Brain)"))
sub$CTs <- "Neuron"
DefaultAssay(sub) <- "SCT"
sub$replicate <- sub$orig.ident
set.seed(1234)
sub$replicate[which(sub$orig.ident == "snuc")] <- sample(paste0("sNuc",1:2),length(which(sub$orig.ident == "snuc")),replace = T)
table(sub$replicate)
for(assay_name in Assays(sub)) {
  sub[[assay_name]] <- as(sub[[assay_name]], "Assay5")
}
excitNeuron.diff <- Libra::run_de(input = sub,
                                  meta = sub@meta.data,
                                  cell_type_col = "CTs",
                                  replicate_col = "replicate",
                                  label_col = "Origin",
                                  de_family = "pseudobulk",
                                  de_method = "edgeR",
                                  de_type = "QLF") #Negative FC means higher in this study

sub <- subset(combined.sct, subset = CTOrigin %in% c("Neuron (This study)","Inhibitory Neuron (Brain)"))
sub$CTs <- "Neuron"
DefaultAssay(sub) <- "SCT"
sub$replicate <- sub$orig.ident
set.seed(1234)
sub$replicate[which(sub$orig.ident == "snuc")] <- sample(paste0("sNuc",1:2),length(which(sub$orig.ident == "snuc")),replace = T)
table(sub$replicate)
for(assay_name in Assays(sub)) {
  sub[[assay_name]] <- as(sub[[assay_name]], "Assay5")
}
inhibNeuron.diff <- Libra::run_de(input = sub,
                                  meta = sub@meta.data,
                                  cell_type_col = "CTs",
                                  replicate_col = "replicate",
                                  label_col = "Origin",
                                  de_family = "pseudobulk",
                                  de_method = "edgeR",
                                  de_type = "QLF") #Negative FC means higher in this study

sub <- subset(combined.sct, subset = CTOrigin %in% c("OPC (This study)","OPC (Brain)"))
sub$CTs <- "OPC"
DefaultAssay(sub) <- "SCT"
sub$replicate <- sub$orig.ident
set.seed(1234)
sub$replicate[which(sub$orig.ident == "snuc")] <- sample(paste0("sNuc",1:2),length(which(sub$orig.ident == "snuc")),replace = T)
table(sub$replicate)
for(assay_name in Assays(sub)) {
  sub[[assay_name]] <- as(sub[[assay_name]], "Assay5")
}
opc.diff <- Libra::run_de(input = sub,
                          meta = sub@meta.data,
                          cell_type_col = "CTs",
                          replicate_col = "replicate",
                          label_col = "Origin",
                          de_family = "pseudobulk",
                          de_method = "edgeR",
                          de_type = "QLF") #Negative FC means higher in this study

save(astrocyte.diff,microglia.diff,inhibNeuron.diff,excitNeuron.diff,opc.diff,file = "./OtherStudies/GSE179590/sNucIntegration/DifferentialExpressionResults.RData",compression_level = 6)
load("./OtherStudies/GSE179590/sNucIntegration/DifferentialExpressionResults.RData")
library(WebGestaltR)
inp_microglia_list <- list("This study" = microglia.diff$gene[microglia.diff$p_val_adj <= 0.01 & microglia.diff$avg_logFC > 0.5],
                           "Brain" = microglia.diff$gene[microglia.diff$p_val_adj <= 0.01 & microglia.diff$avg_logFC < -0.5])

micBrain <- WebGestaltR(enrichMethod = "ORA",
                        organism = "hsapiens",
                        enrichDatabase = "geneontology_Biological_Process",
                        interestGene = inp_microglia_list$Brain,
                        interestGeneType = "genesymbol",
                        referenceGene = rownames(combined.sct@assays$RNA),
                        referenceGeneType = "genesymbol",
                        sigMethod = "fdr",
                        fdrMethod = "BH",
                        fdrThr = 0.01,
                        isOutput = F
)

micThis <- WebGestaltR(enrichMethod = "ORA",
                       organism = "hsapiens",
                       enrichDatabase = "geneontology_Biological_Process",
                       interestGene = inp_microglia_list$`This study`,
                       interestGeneType = "genesymbol",
                       referenceGene = rownames(combined.sct@assays$RNA),
                       referenceGeneType = "genesymbol",
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
                             interestGeneType = "genesymbol",
                             referenceGene = rownames(combined.sct@assays$RNA),
                             referenceGeneType = "genesymbol",
                             sigMethod = "fdr",
                             fdrMethod = "BH",
                             fdrThr = 0.01,
                             isOutput = F
)

excitNeuThis <- WebGestaltR(enrichMethod = "ORA",
                            organism = "hsapiens",
                            enrichDatabase = "geneontology_Biological_Process",
                            interestGene = inp_excitneuron_list$`This study`,
                            interestGeneType = "genesymbol",
                            referenceGene = rownames(combined.sct@assays$RNA),
                            referenceGeneType = "genesymbol",
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
                             interestGeneType = "genesymbol",
                             referenceGene = rownames(combined.sct@assays$RNA),
                             referenceGeneType = "genesymbol",
                             sigMethod = "fdr",
                             fdrMethod = "BH",
                             fdrThr = 0.01,
                             isOutput = F
)

inhibNeuThis <- WebGestaltR(enrichMethod = "ORA",
                            organism = "hsapiens",
                            enrichDatabase = "geneontology_Biological_Process",
                            interestGene = inp_inhibneuron_list$`This study`,
                            interestGeneType = "genesymbol",
                            referenceGene = rownames(combined.sct@assays$RNA),
                            referenceGeneType = "genesymbol",
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
                          interestGeneType = "genesymbol",
                          referenceGene = rownames(combined.sct@assays$RNA),
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
                         referenceGene = rownames(combined.sct@assays$RNA),
                         referenceGeneType = "genesymbol",
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
                        interestGeneType = "genesymbol",
                        referenceGene = rownames(combined.sct@assays$RNA),
                        referenceGeneType = "genesymbol",
                        sigMethod = "fdr",
                        fdrMethod = "BH",
                        fdrThr = 0.01,
                        isOutput = F
)

opcThis <- WebGestaltR(enrichMethod = "ORA",
                       organism = "hsapiens",
                       enrichDatabase = "geneontology_Biological_Process",
                       interestGene = inp_opc_list$`This study`,
                       interestGeneType = "genesymbol",
                       referenceGene = rownames(combined.sct@assays$RNA),
                       referenceGeneType = "genesymbol",
                       sigMethod = "fdr",
                       fdrMethod = "BH",
                       fdrThr = 0.01,
                       isOutput = F
)

save(astroBrain,astroThis,micBrain,micThis,excitNeuBrain,excitNeuThis,inhibNeuBrain,inhibNeuThis,opcThis,opcBrain,file = "./OtherStudies/GSE179590/sNucIntegration/GODifferences.RData")

load("./OtherStudies/GSE179590/sNucIntegration/DifferentialExpressionResults.RData",verbose = T)
write.table(microglia.diff[which(abs(microglia.diff$avg_logFC) > 0.5 & microglia.diff$p_val_adj <= 0.01),],"./OtherStudies/GSE179590/sNucIntegration/DiffExp_Microglia.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(astrocyte.diff[which(abs(astrocyte.diff$avg_logFC) > 0.5 & astrocyte.diff$p_val_adj <= 0.01),],"./OtherStudies/GSE179590/sNucIntegration/DiffExp_Astrocyte.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(inhibNeuron.diff[which(abs(inhibNeuron.diff$avg_logFC) > 0.5 & inhibNeuron.diff$p_val_adj <= 0.01),],"./OtherStudies/GSE179590/sNucIntegration/DiffExp_InhibNeuron.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(excitNeuron.diff[which(abs(excitNeuron.diff$avg_logFC) > 0.5 & excitNeuron.diff$p_val_adj <= 0.01),],"./OtherStudies/GSE179590/sNucIntegration/DiffExp_ExcitNeuron.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(opc.diff[which(abs(opc.diff$avg_logFC) > 0.5 & opc.diff$p_val_adj <= 0.01),],"./OtherStudies/GSE179590/sNucIntegration/DiffExp_OPC.txt",sep = "\t",row.names = F, col.names = T, quote = F)

load("./OtherStudies/GSE179590/sNucIntegration/GODifferences.RData",verbose = T)
write.table(astroBrain,file = "./OtherStudies/GSE179590/sNucIntegration/GO_AstrocyteBrain.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(astroThis,file = "./OtherStudies/GSE179590/sNucIntegration/GO_AstrocyteThis.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(micBrain,file = "./OtherStudies/GSE179590/sNucIntegration/GO_MicrogliaBrain.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(micThis,file = "./OtherStudies/GSE179590/sNucIntegration/GO_MicrogliaThis.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(excitNeuBrain,file = "./OtherStudies/GSE179590/sNucIntegration/GO_ExcitNeuBrain.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(excitNeuThis,file = "./OtherStudies/GSE179590/sNucIntegration/GO_ExcitNeuThis.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(inhibNeuBrain,file = "./OtherStudies/GSE179590/sNucIntegration/GO_InhibNeuBrain.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(inhibNeuThis,file = "./OtherStudies/GSE179590/sNucIntegration/GO_InhibNeuThis.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(opcBrain,file = "./OtherStudies/GSE179590/sNucIntegration/GO_OPCBrain.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(opcThis,file = "./OtherStudies/GSE179590/sNucIntegration/GO_OPCThis.txt",sep = "\t",row.names = F, col.names = T, quote = F)

#Linnarsson
combined.sct <- readRDS("./OtherStudies/Linnarsson/Combined_SeuratObj.rds")
table(combined.sct$CTOrigin)


sub <- subset(combined.sct, subset = CTOrigin %in% c("Microglia (This study)","Microglia (Brain)"))
DefaultAssay(sub) <- "SCT"
sub$replicate <- sub$donor_id
set.seed(1234)
sub$replicate[which(sub$donor_id == "sNuc")] <- sample(paste0("sNuc",1:2),length(which(sub$donor_id == "sNuc")),replace = T)
table(sub$replicate)
for(assay_name in Assays(sub)) {
  sub[[assay_name]] <- as(sub[[assay_name]], "Assay5")
}
sub[["RNA"]] <- JoinLayers(sub[["RNA"]])
microglia.diff <- Libra::run_de(input = sub,
                                meta = sub@meta.data,
                                cell_type_col = "CTs",
                                replicate_col = "replicate",
                                label_col = "Origin",
                                de_family = "pseudobulk",
                                de_method = "edgeR",
                                de_type = "QLF") #Negative FC means higher in this study

sub <- subset(combined.sct, subset = CTOrigin %in% c("Astrocyte (This study)","Astrocyte (Brain)") & donor_id %in% c("H18.30.002", "H19.30.001", "H19.30.002", "sNuc"))
DefaultAssay(sub) <- "SCT"
sub$replicate <- sub$donor_id
set.seed(1234)
sub$replicate[which(sub$donor_id == "sNuc")] <- sample(paste0("sNuc",1:20),length(which(sub$donor_id == "sNuc")),replace = T)
table(sub$replicate)
for(assay_name in Assays(sub)) {
  sub[[assay_name]] <- as(sub[[assay_name]], "Assay5")
}
sub[["RNA"]] <- JoinLayers(sub[["RNA"]])
astrocyte.diff <- Libra::run_de(input = sub,
                                meta = sub@meta.data,
                                cell_type_col = "CTs",
                                replicate_col = "replicate",
                                label_col = "Origin",
                                de_family = "pseudobulk",
                                de_method = "edgeR",
                                de_type = "QLF") #Negative FC means higher in this study

sub <- subset(combined.sct, subset = CTOrigin %in% c("Neuron (This study)","Glutamatergic Neuron (Brain)"))
sub$CTs <- "Neuron"
DefaultAssay(sub) <- "SCT"
sub$replicate <- sub$donor_id
set.seed(1234)
sub$replicate[which(sub$donor_id == "sNuc")] <- sample(paste0("sNuc",1:2),length(which(sub$donor_id == "sNuc")),replace = T)
table(sub$replicate)
for(assay_name in Assays(sub)) {
  sub[[assay_name]] <- as(sub[[assay_name]], "Assay5")
}
sub[["RNA"]] <- JoinLayers(sub[["RNA"]])
excitNeuron.diff <- Libra::run_de(input = sub,
                                  meta = sub@meta.data,
                                  cell_type_col = "CTs",
                                  replicate_col = "replicate",
                                  label_col = "Origin",
                                  de_family = "pseudobulk",
                                  de_method = "edgeR",
                                  de_type = "QLF") #Negative FC means higher in this study

sub <- subset(combined.sct, subset = CTOrigin %in% c("Neuron (This study)","GABAergic Neuron (Brain)"))
sub$CTs <- "Neuron"
DefaultAssay(sub) <- "SCT"
sub$replicate <- sub$donor_id
set.seed(1234)
sub$replicate[which(sub$donor_id == "sNuc")] <- sample(paste0("sNuc",1:2),length(which(sub$donor_id == "sNuc")),replace = T)
table(sub$replicate)
for(assay_name in Assays(sub)) {
  sub[[assay_name]] <- as(sub[[assay_name]], "Assay5")
}
sub[["RNA"]] <- JoinLayers(sub[["RNA"]])
inhibNeuron.diff <- Libra::run_de(input = sub,
                                  meta = sub@meta.data,
                                  cell_type_col = "CTs",
                                  replicate_col = "replicate",
                                  label_col = "Origin",
                                  de_family = "pseudobulk",
                                  de_method = "edgeR",
                                  de_type = "QLF") #Negative FC means higher in this study

sub <- subset(combined.sct, subset = CTOrigin %in% c("OPC (This study)","OPC (Brain)") & donor_id %in% c("H18.30.002", "H19.30.001", "H19.30.002", "sNuc"))
sub$CTs <- "OPC"
DefaultAssay(sub) <- "SCT"
sub$replicate <- sub$donor_id
set.seed(1234)
sub$replicate[which(sub$donor_id == "sNuc")] <- sample(paste0("sNuc",1:2),length(which(sub$donor_id == "sNuc")),replace = T)
table(sub$replicate)
for(assay_name in Assays(sub)) {
  sub[[assay_name]] <- as(sub[[assay_name]], "Assay5")
}
sub[["RNA"]] <- JoinLayers(sub[["RNA"]])
opc.diff <- Libra::run_de(input = sub,
                          meta = sub@meta.data,
                          cell_type_col = "CTs",
                          replicate_col = "replicate",
                          label_col = "Origin",
                          de_family = "pseudobulk",
                          de_method = "edgeR",
                          de_type = "QLF") #Negative FC means higher in this study

save(astrocyte.diff,microglia.diff,inhibNeuron.diff,excitNeuron.diff,opc.diff,file = "./OtherStudies/Linnarsson/DifferentialExpressionResults.RData",compression_level = 6)
load("./OtherStudies/Linnarsson/DifferentialExpressionResults.RData")
library(WebGestaltR)
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

save(astroBrain,astroThis,micBrain,micThis,excitNeuBrain,excitNeuThis,inhibNeuBrain,inhibNeuThis,opcThis,opcBrain,file = "./OtherStudies/Linnarsson/GODifferences.RData")

load("./OtherStudies/Linnarsson/DifferentialExpressionResults.RData",verbose = T)
write.table(microglia.diff[which(abs(microglia.diff$avg_logFC) > 0.5 & microglia.diff$p_val_adj <= 0.01),],"./OtherStudies/Linnarsson/DiffExp_Microglia.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(astrocyte.diff[which(abs(astrocyte.diff$avg_logFC) > 0.5 & astrocyte.diff$p_val_adj <= 0.01),],"./OtherStudies/Linnarsson/DiffExp_Astrocyte.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(inhibNeuron.diff[which(abs(inhibNeuron.diff$avg_logFC) > 0.5 & inhibNeuron.diff$p_val_adj <= 0.01),],"./OtherStudies/Linnarsson/DiffExp_InhibNeuron.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(excitNeuron.diff[which(abs(excitNeuron.diff$avg_logFC) > 0.5 & excitNeuron.diff$p_val_adj <= 0.01),],"./OtherStudies/Linnarsson/DiffExp_ExcitNeuron.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(opc.diff[which(abs(opc.diff$avg_logFC) > 0.5 & opc.diff$p_val_adj <= 0.01),],"./OtherStudies/Linnarsson/DiffExp_OPC.txt",sep = "\t",row.names = F, col.names = T, quote = F)

load("./OtherStudies/Linnarsson/GODifferences.RData",verbose = T)
write.table(astroBrain,file = "./OtherStudies/Linnarsson/GO_AstrocyteBrain.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(astroThis,file = "./OtherStudies/Linnarsson/GO_AstrocyteThis.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(micBrain,file = "./OtherStudies/Linnarsson/GO_MicrogliaBrain.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(micThis,file = "./OtherStudies/Linnarsson/GO_MicrogliaThis.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(excitNeuBrain,file = "./OtherStudies/Linnarsson/GO_ExcitNeuBrain.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(excitNeuThis,file = "./OtherStudies/Linnarsson/GO_ExcitNeuThis.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(inhibNeuBrain,file = "./OtherStudies/Linnarsson/GO_InhibNeuBrain.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(inhibNeuThis,file = "./OtherStudies/Linnarsson/GO_InhibNeuThis.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(opcBrain,file = "./OtherStudies/Linnarsson/GO_OPCBrain.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(opcThis,file = "./OtherStudies/Linnarsson/GO_OPCThis.txt",sep = "\t",row.names = F, col.names = T, quote = F)

