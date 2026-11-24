###############################################################################
# This script generates the differential expression analysis results of
# astrocytes, radial glia and neurons between organoids and the triculture.
# 
# Results are shown in Supplementary Table 4.5 of the manuscript
###############################################################################

library(Seurat)
library(Libra)

load("./IntegrationWithOtherStudies/scOrganoids/Combined_SeuratObj.RData",verbose = T)
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
                                de_type = "QLF")

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
                             de_type = "QLF")

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
                                 de_type = "QLF")

write.table(astrocyte.diff[which(abs(astrocyte.diff$avg_logFC) > 0.5 & astrocyte.diff$p_val_adj <= 0.01),],"./IntegrationWithOtherStudies/scOrganoids/DiffExp_Astrocyte.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(neuron.diff[which(abs(neuron.diff$avg_logFC) > 0.5 & neuron.diff$p_val_adj <= 0.01),],"./IntegrationWithOtherStudies/scOrganoids/DiffExp_Neuron.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(radialglia.diff[which(abs(radialglia.diff$avg_logFC) > 0.5 & radialglia.diff$p_val_adj <= 0.01),],"./IntegrationWithOtherStudies/scOrganoids/DiffExp_RadialGlia.txt",sep = "\t",row.names = F, col.names = T, quote = F)

save(astrocyte.diff,radialglia.diff,neuron.diff,file = "./IntegrationWithOtherStudies/scOrganoids/DiffExp.RData",compression_level = 6)
