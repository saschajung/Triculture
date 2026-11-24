###############################################################################
# This script generates the differential expression analysis results of
# OPCs, astrocytes, radial glia and neurons between the Linnarsson data and the triculture.
# 
# Results are shown in Supplementary Table 4.3 of the manuscript
###############################################################################

library(Seurat)
library(Libra)

combined.sct <- readRDS("./IntegrationWithOtherStudies/Linnarsson/Combined_SeuratObj.rds")

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

save(astrocyte.diff,microglia.diff,inhibNeuron.diff,excitNeuron.diff,opc.diff,file = "./IntegrationWithOtherStudies/Linnarsson/DifferentialExpressionResults.RData",compression_level = 6)

write.table(microglia.diff[which(abs(microglia.diff$avg_logFC) > 0.5 & microglia.diff$p_val_adj <= 0.01),],"./IntegrationWithOtherStudies/Linnarsson/DiffExp_Microglia.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(astrocyte.diff[which(abs(astrocyte.diff$avg_logFC) > 0.5 & astrocyte.diff$p_val_adj <= 0.01),],"./IntegrationWithOtherStudies/Linnarsson/DiffExp_Astrocyte.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(inhibNeuron.diff[which(abs(inhibNeuron.diff$avg_logFC) > 0.5 & inhibNeuron.diff$p_val_adj <= 0.01),],"./IntegrationWithOtherStudies/Linnarsson/DiffExp_InhibNeuron.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(excitNeuron.diff[which(abs(excitNeuron.diff$avg_logFC) > 0.5 & excitNeuron.diff$p_val_adj <= 0.01),],"./IntegrationWithOtherStudies/Linnarsson/DiffExp_ExcitNeuron.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(opc.diff[which(abs(opc.diff$avg_logFC) > 0.5 & opc.diff$p_val_adj <= 0.01),],"./IntegrationWithOtherStudies/Linnarsson/DiffExp_OPC.txt",sep = "\t",row.names = F, col.names = T, quote = F)


