library(CCPlotR)
library(plyr)
library(dplyr)
library(ggraph)
library(ggsci)

source("./aux/colorDefinitions.R")

if(!file.exists("./scRNA-seq/Y1_forPaperFinal/FinalInteractions.txt")){
  source("./aux/CreateInputData_Fig5BDE.R")
}
if(!file.exists("./sNuc-seq/Intercom_forPaper_RNA/FinalInteractions.txt")){
  source("./aux/CreateInputData_Fig5ACE.R")
}

ints_sc <- read.table("./scRNA-seq/Y1_forPaperFinal/FinalInteractions.txt",header = T, sep = "\t")
ints_sn <- read.table("./sNuc-seq/Intercom_forPaper_RNA/FinalInteractions.txt",header = T, sep = "\t")

inp_sc <- data.frame(source = ints_sc$Lig.pop,
                     target = ints_sc$Rec.pop,
                     ligand = ints_sc$Ligand,
                     receptor = ints_sc$Receptor,
                     score = log(ints_sc$score),
                     stringsAsFactors = F)
inp_sc$source[which(inp_sc$source == "Activated Microglia")] <- "Microglia"
inp_sc$target[which(inp_sc$target == "Activated Microglia")] <- "Microglia"
inp_sc$source[which(inp_sc$source == "Proliferating Microglia")] <- "Microglia"
inp_sc$target[which(inp_sc$target == "Proliferating Microglia")] <- "Microglia"
inp_sc <- aggregate(score~.,inp_sc,max)
inp_sc$ID <- paste0(inp_sc$source,"_",inp_sc$target,"_",inp_sc$ligand,"_",inp_sc$receptor)
inp_sc <- inp_sc[which(inp_sc$source %in% c("Microglia","Neuron","Astrocyte") & inp_sc$target %in% c("Microglia","Neuron","Astrocyte")),]

inp_sn <- data.frame(source = ints_sn$Lig.pop,
                     target = ints_sn$Rec.pop,
                     ligand = ints_sn$Ligand,
                     receptor = ints_sn$Receptor,
                     score = log(ints_sn$score),
                     stringsAsFactors = F)
inp_sn$ID <- paste0(inp_sn$source,"_",inp_sn$target,"_",inp_sn$ligand,"_",inp_sn$receptor)
inp_sn <- inp_sn[which(inp_sn$source %in% c("Microglia","Neuron","Astrocyte") & inp_sn$target %in% c("Microglia","Neuron","Astrocyte")),]

commonInts <- intersect(inp_sn$ID,inp_sc$ID)
inp_common <- merge(inp_sc,inp_sn,by = "ID")
inp_common$score <- rowMaxs(as.matrix(inp_common[,c("score.x","score.y")]))
inp_common <- inp_common[,c("source.x","target.x","ligand.x","receptor.x","score")]
colnames(inp_common) <- c("source","target","ligand","receptor","score")

p <- cc_network(inp_common,option = "B",n_top_ints = nrow(inp_common),colours = colorMap,layout = "kk",label_size = 16/.pt, node_size = 2.5) +
  scale_fill_manual(values = colorMap,name = "Cell Type") +
  theme(legend.text=element_text(size=16),
        legend.title = element_text(size = 20))

ggsave("./PaperPlots/Fig5E.png",p,dpi=300,width = 8, height = 8)





