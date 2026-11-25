library(CCPlotR)
library(plyr)
library(dplyr)
library(ggplot2)

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

input_df <- inp_sn %>% mutate(source = factor(source), 
                              target = factor(target)) %>% group_by(source, target, 
                                                                    .drop = FALSE) %>% tally()
brks <- (scales::pretty_breaks(n = 5))(c(round_any(min(input_df$n),10, f = floor), round_any(max(input_df$n), 10, f = ceiling)))
p_sn <- cc_dotplot(inp_sn, option = "A") +
  scale_fill_gradientn(colors = pal_jco_size(5)) +
  scale_size(range = c(10, 20), limits = c(round_any(min(input_df$n),10, f = floor), round_any(max(input_df$n), 10,f = ceiling)), breaks = rev(brks), name = "Number of\ninteractions") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle=45,hjust=1,vjust = 1),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20))
p_sn$guides$guides$size$params$override.aes$fill = rev(pal_jco_size(5))#viridis(5,direction = 1, option = "C")

ggsave("./PaperPlots/Fig5A.png",p_sn,dpi=300,width = 6, height = 6)
