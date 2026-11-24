library(Seurat)

load("./IntegrationWithOtherStudies/GSE179590/Combined_SeuratObj.RData")

obj <- subset(combined.sct,subset = orig.ident != "snuc")
scNuc <- subset(combined.sct,subset = orig.ident == "snuc")

commonGenes <- intersect(rownames(obj@assays$SCT),rownames(scNuc@assays$SCT))

pseudo_ref <- PseudobulkExpression(object = obj,
                                   assay = "SCT",
                                   features = commonGenes,
                                   group.by = "CTs",
                                   method = "aggregate",
                                   layer = "counts")
pseudo_ref$SCT <- apply(pseudo_ref$SCT,2,function(x){(x*1e6)/sum(x)})
cnames <- colnames(pseudo_ref$SCT)
pseudo_ref$SCT <- t(apply(pseudo_ref$SCT,1,function(x){as.numeric(scale(x))}))
colnames(pseudo_ref$SCT) <- cnames

cor(pseudo_ref$SCT,use = "complete.obs")

pseudo_query <- PseudobulkExpression(object = scNuc,
                                     assay = "SCT",
                                     features = commonGenes,
                                     group.by = "CTs",
                                     method = "aggregate",
                                     layer = "counts")
pseudo_query$SCT <- apply(pseudo_query$SCT,2,function(x){(x*1e6)/sum(x)})
cnames <- colnames(pseudo_query$SCT)
pseudo_query$SCT <- t(apply(pseudo_query$SCT,1,function(x){as.numeric(scale(x))}))
colnames(pseudo_query$SCT) <- cnames

cormat <- cor(as.matrix(pseudo_query[["SCT"]]),as.matrix(pseudo_ref[["SCT"]]),method = "pearson",use = "complete.obs")

saveRDS(cormat, file = "./InputData/InputData_Fig4I.rds")
