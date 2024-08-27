########################
# METABOLOMICS FIGURES #
#      HILIC POS       #
########################

# SWIP-10
library(tidyverse)
#install.packages("rjson")
library(rjson)
library(gplots)
library(RColorBrewer)
library(ggbiplot)

# Figures to make
# 1. Heatmap
# 2. PCA
# 3. MWAS
# 4. Pathway analysis

hilicpos <- read_csv("raw_files/swip10_metabolomics/hilicpos/metaboanalyst_hm_200/data_normalized.csv")
names.pos <- hilicpos$X1[-1]

names.mz <- unlist(lapply(str_split(as.character(names.pos), "/"), function(x) as.numeric(x[1])))
names.time <- unlist(lapply(str_split(as.character(names.pos), "/"), function(x) as.numeric(x[2])))

names_rounded <- paste0(round(names.mz,7), "_", round(names.time,7))
names_rounded[1:5]

features <- hilicpos %>% 
    t(.) %>% 
    as.data.frame(.) %>% 
    janitor::row_to_names(1) 

names(features) <- c("Label", names_rounded)
features[1:10,1:5]


########################################
#                MWAS                  #    
########################################

t.test <- read_csv("raw_files/swip10_metabolomics/hilicpos/metaboanalyst_hm_200/t_test.csv") %>% 
    mutate(time = unlist(lapply(str_split(.$X1, "/"), function(x) x[2]))) %>% 
    mutate(time = as.numeric(time)) %>% 
    mutate(posneg = case_when(t.stat >= 0 & p.value < 0.05 ~ 2,
                              t.stat < 0 & p.value < 0.05 ~ 3, #higher in swip10
                              p.value >= 0.05 ~ 1))


png("swip-10/figures/hilicpos_mwas_swip10_p005.png", width =4.5, height = 4, units = 'in', res = 300)

ggplot(data = t.test, aes(x = time, y = `-log10(p)`, color = factor(posneg))) +
    geom_jitter(size = 0.75) +
    #ggtitle("", subtitle = "Features in red are higher in *dhs-11*") +
    scale_color_manual(values=c("grey90", "grey2", "indianred4"), label = c("ns", "higher in N2", "higher in swip-10")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 2)) +
    geom_hline(yintercept = -log10(0.05), color = "black", linetype = "dashed") + #p = 0.05
    #geom_hline(yintercept = -log10(0.10), color = "blue", linetype = "dashed") + #p = 0.1
    theme_classic() +
    ylab("-log10(p)") +
    xlab("retention time") +
    labs(color = "") +
    theme(legend.position = "bottom", plot.subtitle = ggtext::element_markdown()) 

dev.off()

########################################
#                PCA                   #    
########################################

pca.scores <- fromJSON(file = "raw_files/swip10_metabolomics/hilicpos/metaboanalyst_hm_200/pca_score3d_0_.json")
pca.scores

# PCA with own centering and scaling
feat.pca <- features %>% 
    select(contains("_")) %>% 
    map_df(~as.numeric(as.character(.x))) %>% 
    prcomp(., center = T, scale = T)

features$strain <- case_when(features$Label == "swip10" ~ "swip-10",
                             features$Label == "N2" ~ "N2")

png("swip-10/figures/pca_hilicpos.png", width = 4.5, height = 4.5, units = 'in', res = 300)
ggbiplot(feat.pca,ellipse=TRUE,  groups=factor(features$strain), var.axes=F, obs.scale = 1, var.scale = 1)+
    scale_colour_manual(name="Strain", values= c("indianred4", "black"))+
    theme_bw() +
    theme(legend.position = "none") +
    annotate("text", x = 15, y = 44, label = "swip-10", size = 4, fontface = "italic", color = "indianred4") +
    annotate("text", x = -15, y = -30, label = "N2", size = 4, color = "black")
dev.off()  

features[1:5,1:5]


########################################
#              HEATMAP                 #    
########################################
# Heatmap, top 200 features

top200 <- t.test %>%
    arrange(p.value) %>%
    slice(1:200) %>% 
    mutate(mz = unlist(lapply(str_split(.$X1, "/"), function(x) as.numeric(x[1])))) %>% 
    mutate(mz_time = paste0(round(mz,7), "_", round(time,7)))

heat_200 <- features %>% 
    select(top200$mz_time) %>% 
    map_df(~as.character(.x)) %>% 
    map_df(~as.numeric(.x)) 

distCor <- function(x) dist(x, method = "euclidean")
hclustAvg <- function(x) hclust(x, method="ward.D2")

heat <- matrix(data = heat_200 %>% unlist(), nrow = 10) 
labs <- unlist(features$Label) 
labs.swip10 <-  c(rep("N2", 5), rep("swip-10", 5))
row.names(heat) <- labs.swip10
my.pal <- brewer.pal(6, "RdYlBu")

png("swip-10/figures/heatmap_swip10_hilicpos_redyellowblue.png", width = 9, height = 7, units = 'in', res = 300)

heatmap.2(heat,
          Rowv=TRUE,
          Colv=TRUE,
          trace='none',        # turns off trace lines inside the heat map
          density.info="none",
          dendrogram = "none",
          distfun=distCor, 
          hclustfun=hclustAvg,
          #RowSideColors = colSide, 
          #margins = c(5,10),
          scale = "col",
          labRow = c(rep("N2", 5), as.expression(lapply(rownames(heat)[6:10], function(a) bquote(italic(.(a)))))),
          cexRow = 0.8, cexCol = 0.01,
          col = my.pal)     

dev.off()

################
# PATHWAY ANALYSIS
################

# CALL IN PATHWAY FILE
path <- read_tsv("raw_files/swip10_metabolomics/hilicpos/1534265227.46.human_swip10_hilicpos_10000/tables/mcg_pathwayanalysis_human_swip10.tsv")

png("swip-10/figures/hilpos_pathways_bubble_010.png", res= 300, units = "in", width = 4.5, height = 4.5)

path %>% 
    mutate(logp = -log10(`p-value`)) %>% 
    filter(overlap_size > 1) %>% 
    mutate(enrich = overlap_size/pathway_size) %>% 
    mutate(label = paste0(pathway, " (", overlap_size, "/", pathway_size, ")")) %>% 
    filter(`p-value` < 0.1) %>% 
    ggplot(aes(x=logp, y = reorder(label, logp), color = enrich)) +
    geom_point( size = 3) +
    scale_color_gradient(low="grey76", high="grey20")+
    #theme_minimal() +
    xlab("-log10(p-value)") +
    ylab("") +
    theme_bw() +
    #ggtitle("Pathways altered in cat-1 worms", 
    #        subtitle = "Size of bubble represents number of significant hits \nEnrichment is calculated as (Total Hits/Pathway size)") +
    theme(plot.title = element_text(size = 9, face = "bold"),
          plot.subtitle = element_text(size = 7),
          axis.text.x =element_text(size=5),
          axis.text.y =element_text(size=8),
          axis.title=element_text(size=9,face="bold"),
          strip.text = element_text(size=7),
          legend.text=element_text(size=7),
          legend.title=element_text(size=8),
          legend.position=c(0.7,0.18),
          panel.grid.minor = element_blank(),
          legend.key.size = unit(0.5, "lines")) +
    guides(size=guide_legend("Overlap size"))  +
    labs(col = "Enrichment") 
dev.off()

