library(ggpubr)
library(colorspace)

if(!file.exists("./InputData/InputData_Fig4H.rds")){
  source("./aux/CreateInputData_Fig4H.R")
}

neuronThis_SynGO <- readRDS("./InputData/InputData_Fig4H.rds")
neuronThis_SynGO$description[which(neuronThis_SynGO$description == "modulation of chemical synaptic transmission")] <- "modulation of chemical\nsynaptic transmission" 
neuronThis_SynGO$description[which(neuronThis_SynGO$description == "regulation of synapse organization")] <- "regulation of\nsynapse organization" 

p <- ggbarplot(neuronThis_SynGO,
               x = "description",
               y = "enrichmentRatio",
               fill = "FDR",
               sort.val = "asc") +
  guides(fill = guide_colorbar(barheight = 1, barwidth = 15)) +
  scale_fill_gradient(high = lighten("#EFC000FF", amount = 0.8),low="#EFC000FF")+
  coord_flip() +
  xlab("") + 
  ylab("Enrichment Ratio") +
  theme(legend.title = element_text(size = 24),
        legend.text = element_text(size = 22),
        axis.text = element_text(size = 22),
        axis.title = element_text(size = 24))

ggsave("./PaperPlots/Fig4H.png",p,dpi = 300,width = 8, height = 8)
