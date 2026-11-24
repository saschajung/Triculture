library(ggheatmap)
library(ggtree)

source("./aux/colorDefinitions.R")

if(!file.exists("./InputData/InputData_Fig4J.rds")){
  source("./aux/CreateInputData_Fig4J.R")
}

cormat <- readRDS("./InputData/InputData_Fig4J.rds")
p <- ggheatmap(cormat,
               cluster_rows = T,
               cluster_cols = T,
               hclust_method = "ward.D2",
               dist_method = "maximum",
               border = "black")
p$plotlist[[1]] <- p$plotlist[[1]] + 
  geom_tile(color = "black",linewidth = 1) +
  scale_fill_gradientn(colors = pal_jco_heatmap(3),values = scales::rescale(c(min(cormat),0,max(cormat)),to = c(0,1) )) +
  theme(axis.text = element_text(size = 16, color = "black"),
        legend.text = element_text(size = 16, hjust = 1),
        legend.title = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
  ) +
  labs(fill = "Correlation")
p$plotlist[[2]] <- p$plotlist[[2]] +
  geom_tree(linewidth = 1)
p$plotlist[[3]] <- p$plotlist[[3]] +
  geom_tree(linewidth = 1)
p

ggsave(filename = "./PaperPlots/Fig4J.png",plot = p, dpi = 300, width = 8, height = 6)
