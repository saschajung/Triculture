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

input_df <- inp_sc %>% mutate(source = factor(source), 
                              target = factor(target)) %>% group_by(source, target, 
                                                                    .drop = FALSE) %>% tally()
brks <- (scales::pretty_breaks(n = 5))(c(round_any(min(input_df$n),10, f = floor), round_any(max(input_df$n), 10, f = ceiling)))
p_sc <- cc_dotplot(inp_sc, option = "A") +
  scale_fill_gradientn(colors = pal_jco_size(6)) +
  scale_size(range = c(10, 20), limits = c(round_any(min(input_df$n),10, f = floor), round_any(max(input_df$n), 10,f = ceiling)), breaks = rev(brks), name = "Number of\ninteractions") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle=45,hjust=1,vjust = 1),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20))
p_sc$guides$guides$size$params$override.aes$fill = rev(pal_jco_size(6))#viridis(5,direction = 1, option = "C")

ggsave("./PaperPlots/Fig5B.png",p_sc,dpi=300,width = 6, height = 6)
