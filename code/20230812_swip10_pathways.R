##################
# Pathway analysis with both columns
##################
library(tidyverse)

hilic.path <- read_tsv("raw_files/swip10_metabolomics/hilicpos/1534265227.46.human_swip10_hilicpos_10000/tables/mcg_pathwayanalysis_human_swip10.tsv") %>% 
    mutate(column = "HILIC (positive)") %>% 
    filter(`p-value` < 0.1)

c18.path <- read_tsv("raw_files/swip10_metabolomics/c18neg/1604939205.08.swip10_human_correct/tables/mcg_pathwayanalysis_swip10_human_correct.tsv") %>% 
    mutate(column = "C18 (negative)") %>% 
    filter(`p-value` < 0.1)

png("swip-10/figures/metabolomics/hilicpos/20231122_hilic_pathways_010.png", res= 300, units = "in", width = 4.8, height = 4.5) #Save image to disc as png

hilic.path %>% 
    mutate(logp = -log10(`p-value`)) %>% 
    mutate(enrich = overlap_size/pathway_size) %>% 
    filter(overlap_size > 1) %>% 
    mutate(label = paste0(pathway, " (", overlap_size, "/", pathway_size, ")")) %>% 
    ggplot(aes(x=logp, y = reorder(label, logp), color = enrich)) +
    geom_point( size = 3) +
    scale_color_gradient(low="royalblue3", high="maroon")+ # Color of bubbles. Will take hex codes
    xlab("-log10(p-value)") +
    ylab("") +
    theme_bw() +
    theme(plot.title = element_text(size = 9, face = "bold"),
          plot.subtitle = element_text(size = 7),
          axis.text.x =element_text(size=8, color= "black"), # Change size of x-axis labels
          axis.text.y =element_text(size=9, color = "black"), # Change size of y-axis labels
          axis.title=element_text(size=9,face="bold"),
          strip.text = element_text(size=10),
          legend.text=element_text(size=7),
          legend.title=element_text(size=8),
          legend.position=c(0.65,0.14),
          panel.grid.minor = element_blank(),
          legend.key.size = unit(0.5, "lines")) +
    guides(size=guide_legend("Overlap size"))  +
    labs(col = "Enrichment") +
    facet_wrap(~column, nrow = 1)

dev.off()

png("swip-10/figures/metabolomics/c18neg/20231122_c18_pathways_010.png", res= 300, units = "in", width = 4.8, height = 4.5) #Save image to disc as png

c18.path %>% 
    mutate(logp = -log10(`p-value`)) %>% 
    mutate(enrich = overlap_size/pathway_size) %>% 
    filter(overlap_size > 1) %>% 
    mutate(label = paste0(pathway, " (", overlap_size, "/", pathway_size, ")")) %>% 
    ggplot(aes(x=logp, y = reorder(label, logp), color = enrich)) +
    geom_point( size = 3) +
    scale_color_gradient(low="royalblue3", high="maroon")+ # Color of bubbles. Will take hex codes
    xlab("-log10(p-value)") +
    ylab("") +
    theme_bw() +
    theme(plot.title = element_text(size = 9, face = "bold"),
          plot.subtitle = element_text(size = 7),
          axis.text.x =element_text(size=8, color= "black"), # Change size of x-axis labels
          axis.text.y =element_text(size=9, color = "black"), # Change size of y-axis labels
          axis.title=element_text(size=9,face="bold"),
          strip.text = element_text(size=10),
          legend.text=element_text(size=7),
          legend.title=element_text(size=8),
          legend.position=c(0.65,0.14),
          panel.grid.minor = element_blank(),
          legend.key.size = unit(0.5, "lines")) +
    guides(size=guide_legend("Overlap size"))  +
    labs(col = "Enrichment") +
    facet_wrap(~column, nrow = 1)

dev.off()

