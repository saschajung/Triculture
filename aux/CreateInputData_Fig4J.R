library(Seurat)

load("./IntegrationWithOtherStudies/scOrganoids/Combined_SeuratObj.RData",verbose = T)

sub <- subset(combined.sct, subset = orig.ident != "Y1")
y1 <- subset(combined.sct, subset = orig.ident == "Y1")

rm(combined.sct)

commonGenes <- intersect(rownames(sub@assays$SCT),rownames(y1@assays$SCT))

pseudo_ref <- PseudobulkExpression(object = sub,
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

pseudo_query <- PseudobulkExpression(object = y1,
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

saveRDS(cormat, file = "./InputData/InputData_Fig4J.rds")

rm(y1)
rm(sub)

