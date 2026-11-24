library(ggpubr)
library(CCPlotR)
library(plyr)
library(dplyr)

source("./aux/colorDefinitions.R")

if(!file.exists("./scRNA-seq/Y1_forPaperFinal/FinalInteractions.txt")){
  source("./aux/CreateInputData_Fig5BDE.R")
}

ints_sc <- read.table("./scRNA-seq/Y1_forPaperFinal/FinalInteractions.txt",header = T, sep = "\t")

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

df <- as.data.frame(table(inp_sc$target,inp_sc$receptor))
df <- df[which(df$Freq > 0),]
df$Var1 <- factor(as.character(df$Var1),levels = c("Neuron","Astrocyte","Microglia"))
df$Var2 <- factor(as.character(df$Var2), levels = unique(as.character(df$Var2)))
df$Var2_unique <- paste0(df$Var2, "_", seq_len(nrow(df)))
rownames(df) <- NULL

oneCT_sc <- names(which(table(df$Var2)==1))
multInts_sc <- names(which(table(inp_sc$receptor) > 1))
oneInt_sc <- names(which(table(inp_sc$receptor) == 1))
oneCTmultInts_sc <- intersect(oneCT_sc,multInts_sc)
oneCToneInt_sc <- intersect(oneCT_sc,oneInt_sc)

p <- ggbarplot(df, x = "Var2_unique", y = "Freq",
               fill = "Var1",               # change fill color by cyl
               color = "white",            # Set bar border colors to white
               palette = "jco",            # jco journal color palett. see ?ggpar
               sort.val = "desc",          # Sort the value in dscending order
               sort.by.groups = T,     # Don't sort inside each group
               x.text.angle = 90,           # Rotate vertically x axis texts
               xlab = "Receptor"
) +
  font("legend.title",size = 20) +
  font("legend.text",size = 16) +
  font("xy.title", size = 20) +
  font("xy.text",size = 16) +
  scale_fill_manual(values = colorMap,name = "Cell Type")
pb <- ggplot_build(p)

# Access x-axis data (after sorting)
ordered_levels <- pb$data[[1]]$x
ordered_labels <- pb$layout$panel_params[[1]]$x$get_labels()
new_labels <- gsub("_[0-9]+","",ordered_labels)

p <- p + scale_x_discrete(labels = new_labels)
p

ggsave("./PaperPlots/Fig5D.png",p,dpi=300,width = 8, height = 6)
