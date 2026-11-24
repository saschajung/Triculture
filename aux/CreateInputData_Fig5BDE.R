library(InterCom)

scNuc <- readRDS("./sNuc-seq/sNuc_SeuratObj.rds")

inpmat <- as.matrix(scNuc@assays$RNA@counts)
anno.tbl <- data.frame(cell.name = rownames(scNuc@meta.data),
                       cell.type = unname(scNuc$CTs),
                       stringsAsFactors = F)
anno.tbl$cell.type <- gsub(" ","\\.",anno.tbl$cell.type)
anno.tbl$cell.type <- gsub("-","\\.",anno.tbl$cell.type)
anno.tbl$cell.type <- gsub("\\+","\\.",anno.tbl$cell.type)
anno.tbl$cell.name <- gsub("-","\\.",anno.tbl$cell.name)
colnames(inpmat) <- gsub("-","\\.",colnames(inpmat))
inpmat[1:5,1:5]
InterCom(data=inpmat,
         anno.tbl = anno.tbl,
         species = "HUMAN",
         ncores = 4,
         sig.cutoff = 0,
         z.score.cutoff = 0,
         tissue.name = "forPaper_RNA",
         out.path = "./sNuc-seq/Intercom_forPaper_RNA/")


load("./sNuc-seq/Intercom_forPaper_RNA/output_forPaper_RNA.RData",verbose = T)
ints <- output$final
ints <- ints[ints$score >= 0.9,]
pop_map <- c("OLIG1.GFAP..Cells" = "OLIG1+GFAP+ Cells",
             "Intermediate.Astrocyte.Progenitor" = "Intermediate Astrocyte Progenitor",
             "Radial.Glia" = "Radial Glia",
             "Fibroblast.like" = "Fibroblast-like",
             "Astrocyte" = "Astrocyte",
             "Neuron" = "Neuron",
             "OPC" = "OPC",
             "Microglia" = "Microglia")
ints$Lig.pop <- pop_map[ints$Lig.pop]
ints$Rec.pop <- pop_map[ints$Rec.pop]

write.table(ints,file = "./sNuc-seq/Intercom_forPaper_RNA/FinalInteractions.txt",sep = "\t",row.names = F,col.names = T,quote = F)
