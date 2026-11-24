library(ggpubr)

source("./aux/colorDefinitions.R")

if(!file.exists("./sNuc-seq/Intercom_forPaper_RNA/FinalInteractions.txt")){
  source("./aux/CreateInputData_Fig5ACE.R")
}

ints_sn <- read.table("./sNuc-seq/Intercom_forPaper_RNA/FinalInteractions.txt",header = T, sep = "\t")

inp_sn <- data.frame(source = ints_sn$Lig.pop,
                     target = ints_sn$Rec.pop,
                     ligand = ints_sn$Ligand,
                     receptor = ints_sn$Receptor,
                     score = log(ints_sn$score),
                     stringsAsFactors = F)
inp_sn$ID <- paste0(inp_sn$source,"_",inp_sn$target,"_",inp_sn$ligand,"_",inp_sn$receptor)
inp_sn <- inp_sn[which(inp_sn$source %in% c("Microglia","Neuron","Astrocyte") & inp_sn$target %in% c("Microglia","Neuron","Astrocyte")),]

df <- as.data.frame(table(inp_sn$target,inp_sn$receptor))
df <- df[which(df$Freq > 0),]
df$Var1 <- factor(as.character(df$Var1),levels = c("Neuron","Astrocyte","Microglia"))
df$Var2 <- factor(as.character(df$Var2), levels = unique(as.character(df$Var2)))
df$Var2_unique <- paste0(df$Var2, "_", seq_len(nrow(df)))
rownames(df) <- NULL

oneCT_sn <- names(which(table(df$Var2)==1))
multInts_sn <- names(which(table(inp_sn$receptor) > 1))
oneInt_sn <- names(which(table(inp_sn$receptor) == 1))
oneCTmultInts_sn <- intersect(oneCT_sn,multInts_sn)
oneCToneInt_sn <- intersect(oneCT_sn,oneInt_sn)

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

ggsave("./PaperPlots/Fig5C.png",p,dpi=300,width = 8, height = 6)


