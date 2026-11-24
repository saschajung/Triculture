###############################################################################
# This script generates the differential expression analysis results of
# OPCs, astrocytes, radial glia and neurons between the data in GSE179590 and the triculture.
# 
# Results are shown in Supplementary Table 4.1 of the manuscript
###############################################################################

library(Seurat)
library(Libra)

load("./IntegrationWithOtherStudies/GSE179590/Combined_SeuratObj.RData",verbose = T)

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
                                de_type = "QLF")

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
                                de_type = "QLF")

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
                                  de_type = "QLF")

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
                                  de_type = "QLF")

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
                          de_type = "QLF")

save(astrocyte.diff,microglia.diff,inhibNeuron.diff,excitNeuron.diff,opc.diff,file = "./IntegrationWithOtherStudies/GSE179590/DifferentialExpressionResults.RData",compression_level = 6)

write.table(microglia.diff[which(abs(microglia.diff$avg_logFC) > 0.5 & microglia.diff$p_val_adj <= 0.01),],"./IntegrationWithOtherStudies/GSE179590/DiffExp_Microglia.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(astrocyte.diff[which(abs(astrocyte.diff$avg_logFC) > 0.5 & astrocyte.diff$p_val_adj <= 0.01),],"./IntegrationWithOtherStudies/GSE179590/DiffExp_Astrocyte.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(inhibNeuron.diff[which(abs(inhibNeuron.diff$avg_logFC) > 0.5 & inhibNeuron.diff$p_val_adj <= 0.01),],"./IntegrationWithOtherStudies/GSE179590/DiffExp_InhibNeuron.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(excitNeuron.diff[which(abs(excitNeuron.diff$avg_logFC) > 0.5 & excitNeuron.diff$p_val_adj <= 0.01),],"./IntegrationWithOtherStudies/GSE179590/DiffExp_ExcitNeuron.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(opc.diff[which(abs(opc.diff$avg_logFC) > 0.5 & opc.diff$p_val_adj <= 0.01),],"./IntegrationWithOtherStudies/GSE179590/DiffExp_OPC.txt",sep = "\t",row.names = F, col.names = T, quote = F)

