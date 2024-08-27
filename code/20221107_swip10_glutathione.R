####
# Determining levels of glutathione in SWIP mutants, merged from two runs
# 05.13.2020
# Units of concentration
## cyss, cys, gsh, gssg: uM
## gsh_prot, gssg_prot: nmol/mg total protein
## gsh_intracell, gssg_intracell: uM
## eh_ggs/gsh: ratio, unitless
## protein: ug/ul
####

library(tidyverse)
library(cowplot)

dat <- read.csv("raw_files/glutathione/gsh_swips_october_may_merged.csv", header = T, sep =",")

dat.plot.swip10 <- dat %>% 
    mutate(ratio = gsh_prot/gssg_prot) %>% 
    filter(strain != "swip20") %>% 
    mutate(strain.plot = case_when(strain == "swip10" ~ "swip10",
                                   strain == "N2" ~ "N2")) %>% 
    mutate(strain.plot = fct_reorder(strain.plot, eh_gssg.gsh, .desc = F))

# Statistical tests
swip10.val <- dat.plot.swip10 %>% 
    filter(strain == "swip10")
n2.val <- dat.plot.swip10 %>% 
    filter(strain == "N2")

t.eh <- t.test(swip10.val$eh_gssg.gsh, n2.val$eh_gssg.gsh)
t.gsh <- t.test(swip10.val$gsh_prot, n2.val$gsh_prot)
t.gssg <- t.test(swip10.val$gssg_prot, n2.val$gssg_prot)
t.eh # p = 0.019
t.gsh # p = 0.013
t.gssg # p = 0.225

# PLOTS
# Redox potential

#png("swip-10/figures/EHpotential.png", width = 3, height = 4, units = "in", res = 300)
potential <- ggplot(dat.plot.swip10, aes(strain.plot, eh_gssg.gsh, fill = strain.plot)) + 
    geom_boxplot(aes(fill = strain.plot, color = strain.plot), width = 0.5) +
    geom_jitter(aes(color = strain.plot), shape  = 19, width = 0.1) + 
    theme_classic() +
    theme(panel.grid.major = element_blank(), panel.grid.minor =
                      element_blank()) +
    labs(x = "",
         y = "Eh(GSSG/GSH)") +
    scale_fill_manual(values = c("black", "#8B0000")) +
    scale_color_manual(values = c("grey10", "black")) +
    guides(fill="none",
           color = "none") +
    scale_y_continuous(breaks = seq(-270,-80, 20)) +
    theme(plot.title = element_text(size = 10)) +
    scale_x_discrete(labels = c("N2", "*swip-10*")) +
    theme(axis.text.x = ggtext::element_markdown()) +
    geom_segment(aes(x=1, xend=2, y=-80, yend=-80)) +
    geom_segment(aes(x=1, xend=1, y=-80, yend=-83)) +
    geom_segment(aes(x=2, xend=2, y=-80, yend=-83)) +
    annotate("text", x = 1.5, y = -78, label = "*", size = 6)
#potential
#dev.off()

# Reduced glutathione levels
#png("swip-10/figures/gsh.png", width = 3, height = 4, units = "in", res = 300)
gsh <- ggplot(dat.plot.swip10, aes(strain.plot, gsh_prot, fill = strain.plot)) + 
    geom_boxplot(aes(fill = strain.plot, color = strain.plot), width = 0.5) + 
    geom_jitter(aes(color = strain.plot), shape = 19, width = 0.1) +  
    theme_classic() +
    theme(panel.grid.major = element_blank(), panel.grid.minor =
                      element_blank()) +
    labs(x= "",
         y=" Glutathione (GSH) \nnmol/mg total protein") +
    scale_fill_manual(values = c("black", "#8B0000")) +
    scale_color_manual(values = c("grey10", "black")) +
    guides(fill="none",
           color = "none") +
    scale_y_continuous(breaks = seq(0, 40, 10)) +
    theme(plot.title = element_text(size = 10)) +
    scale_x_discrete(labels = c("N2", "*swip-10*")) +
    theme(axis.text.x = ggtext::element_markdown()) +
    geom_segment(aes(x=1, xend=2, y=44, yend=44)) +
    geom_segment(aes(x=1, xend=1, y=44.05, yend=43.3)) +
    geom_segment(aes(x=2, xend=2, y=44.05, yend=43.3)) +
    annotate("text", x = 1.5, y = 45, label = "*", size = 6)
#gsh
#dev.off()
# Oxidized glutathione
#png("swip-10/figures/gssg.png", width = 3, height = 4, units = "in", res = 300)
gssg <- ggplot(dat.plot.swip10, aes(strain.plot, gssg_prot, fill = strain.plot)) + 
    geom_boxplot(aes(fill = strain.plot, color = strain.plot), width = 0.5) + 
    geom_jitter(aes(color = strain.plot), shape = 19, width = 0.1) +  
    theme_classic() +
    theme(panel.grid.major = element_blank(), panel.grid.minor =
              element_blank()) +
    labs(x= "",
         y="Oxidized glutathione (GSSG) \nnmol/mg total protein") +
    scale_fill_manual(values = c("black", "#8B0000")) +
    scale_color_manual(values = c("grey10", "black")) +
    guides(fill="none",
           color = "none") +
    scale_y_continuous(breaks = seq(0, 12, 2)) +
    theme(plot.title = element_text(size = 10)) +
    scale_x_discrete(labels = c("N2", "*swip-10*")) +
    theme(axis.text.x = ggtext::element_markdown())
#gssg
#dev.off()

png("swip-10/figures/gsh_all_arranged.png", width = 6, height = 3, res = 300, units = "in")
plot_grid(potential, gsh, gssg, labels = c('A', 'B', 'C'), label_size = 12, nrow = 1)
dev.off()

dat.plot.swip10 %>% 
    select(strain, gsh_prot, gssg_prot, eh_gssg.gsh) %>% 
    write_csv("swip-10/data/swip10_glutathione.csv")
