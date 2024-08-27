###########
# SEAHORSE DATA: SWIP-10
# black: N2, indianred4: swip-10
###########

library(tidyverse)
#remotes::install_github("wilkelab/ggtext")
library(ggtext)
library(cowplot)

###############
#Seahorse, Koopman et at Protocol
# Basal resp: 5 measurements
# FCCP: 9 measurements
# NaN3: 4 measurements
# 05.03.2018
# Characterizing SWIP mutants: swip10, swip20, and BYN2s
# To access ggplot2, use R 3.2.2
###############

## Creating vector with number of worms to be added to data table
## background well = NA

worms <- c(NA, 30, 62, 90, 100, 130, 
           34,30,85, NA, 100, 120, 
           45, 42, NA, 100, 90, 95, 
           25, 35, 120, 110, 100, NA)

label <- c("BG","swip10", "N2", "N2", "swip20", "swip20", 
           "swip10", "swip10", "N2", "BG", "swip20", "swip20",
           "swip10", "swip10", "BG", "N2", "swip20", "swip20",
           "swip10", "swip10", "N2", "N2", "swip20", "BG")

# length(worms) #repeating worm number to match number of measurements (18)
worm_num <- c(rep(worms[1:24], each = 18))
label_all <- c(rep(label[1:24], each = 18)) 

## Importing data
ocr_raw <- read.csv("raw_files/seahorse/swips_ocr.csv", header = T, sep = ",")

## Adding worm number to the dataframe
ocr_raw$worm_num <- worm_num

## Adding labels
ocr_raw$label <- label_all

## Noramilze to worm number
ocr_raw$norm <- ocr_raw$OCR/ocr_raw$worm_num

## Rearrange data to plot
ocr_worms_plot <- ocr_raw %>% 
    filter(label != "BG") %>% 
    group_by(label, Measurement)%>%
    dplyr::summarise(mean.ocr.p.worm = mean(norm),  sd.ocr.p.worm = sd(norm), count = n()) %>%
    ungroup() %>%
    group_by() %>%
    mutate(group.plot = case_when(label == "swip10" ~ "swip-10",
                                  label == "swip20" ~ "dhs-11",
                                  label == "N2" ~ "N2"))%>%
    mutate(group.plot = fct_reorder(group.plot, mean.ocr.p.worm, .desc = T))

# Save OCR profile plot
#png("swip-10/figures/seahorse_profile.png", width = 5, height = 5, units = "in", res = 300)

profile <- ocr_worms_plot %>% 
    filter(label != "swip20") %>% 
    ggplot(mapping = aes(x = Measurement, y = mean.ocr.p.worm)) + 
    geom_line(aes(color = group.plot), size = 1.04)+
    geom_point(aes(color = group.plot), size = 2.5) +
    ylim(0, 20) +
    labs(x = "Measurement no. (loop)",
         y = "OCR \n(pmol/min per worm)",
         color = "Strain and \ntreatment") +
    scale_color_manual(values = c("black", "#8B0000")) + 
    scale_x_continuous(breaks = seq(1,18,1)) +
    annotate("text",x = 5, y = 15.8, label = "FCCP", size = 4) +
    annotate("text", x = 15, y = 15.8, label = "Sodium azide", size = 4)+
    geom_segment(aes(x=5, xend=5, y=15, yend=13), 
                 arrow = arrow(length = unit(0.25, "cm")), color = "black")+
    geom_segment(aes(x=15, xend=15, y=15, yend= 13), 
                 arrow = arrow(length = unit(0.25, "cm")), color = "black")+
    annotate("text", x = 3.5, y = 4.2, label = 'swip-10', color = "#8B0000", size = 4.5, fontface = "italic") +
    annotate("text", x = 3.5, y = 11, label = "N2", color = "black", size = 4.5) +
    theme_classic() +
    theme(axis.title = element_text(size = 13),
          axis.text = element_text(size = 11),
          legend.position = "none") 
#profile
#dev.off()

### Calculations
# For well labels:
label.well <- rep(c("swip10", "N2", "N2", "swip20", "swip20", 
                    "swip10", "swip10", "N2", "swip20", "swip20",
                    "swip10", "swip10", "N2", "swip20", "swip20",
                    "swip10", "swip10", "N2", "N2", "swip20"), times = 2)


ocr_worms_cal <- ocr_raw %>% 
    filter(label != "BG") %>% 
    filter(Measurement %in% c(1,2,3,8,9,10)) %>%
    mutate(type = case_when(Measurement >=1 & Measurement < 6 ~ "basal",
                            Measurement >=8 & Measurement < 11 ~  "maximal"))  %>% 
    group_by(type, Well)  %>%
    summarise(mean.ocr  =  mean(norm)) %>%
    ungroup() %>%
    #mutate(well.group  = str_sub(Well, 1,1)) %>%
    mutate(label = label.well) %>% 
    mutate(group = case_when(label == "swip10" ~ "swip-10",
                             label == "swip20" ~ "dhs-11",
                             label == "N2" ~ "N2")) %>% 
    mutate(group = fct_reorder(group, mean.ocr, .desc = T)) 

respiratory.measures <- ocr_worms_cal %>% 
    spread(key = type, value = mean.ocr) %>%
    mutate(spare =  maximal - basal)

resp.swip10 <- respiratory.measures %>% 
    filter(label == "swip10")

resp.n2 <- respiratory.measures %>% 
    filter(label == "N2")

t.basal <- t.test(resp.n2$basal, resp.swip10$basal)
t.max <- t.test(resp.n2$maximal, resp.swip10$maximal)
t.spare <- t.test(resp.n2$spare, resp.swip10$spare)

# BASAL RESPIRATION
#png("swip-10/figures/basal.png", width = 3, height = 4, units = "in", res = 300)
basal <- respiratory.measures %>% 
    filter(label != "swip20") %>% 
    ggplot(mapping  = aes(x= group, y  = basal))  +
    geom_boxplot(aes(fill  = group, color = group), width = 0.5)  +
    geom_jitter(aes(color = group), shape  = 19, width = 0.1) +
    scale_fill_manual(values = c("black", "#8B0000")) +
    scale_color_manual(values = c("grey10", "black")) +
    labs(title = "",
         x = "",
         y = "Basal respiration \n(pmol/min per worm)") +
    theme_classic() +
    scale_x_discrete(labels = c("N2", "*swip-10*")) +
    theme(axis.text.x = ggtext::element_markdown(),
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 11),
          legend.position = "none") +
    geom_segment(aes(x=1, xend=2, y=15.05, yend=15.05)) +
    geom_segment(aes(x=1, xend=1, y=15.065, yend=14.85)) +
    geom_segment(aes(x=2, xend=2, y=15.065, yend=14.85)) +
    annotate("text", x = 1.5, y = 15.5, label = "*", size = 6)

#basal
#dev.off()

# MAXIMAL RESPIRATION
#png("swip-10/figures/maximal.png", width = 3, height = 4, units = "in", res = 300)
maximal <- respiratory.measures %>% 
    filter(label != "swip20") %>% 
    ggplot(mapping  = aes(x= group, y  = maximal))  +
    geom_boxplot(aes(fill  = group, color = group), width = 0.5)  +
    geom_jitter(aes(color = group), shape = 19, width = 0.1) +
    scale_fill_manual(values = c("black", "#8B0000")) +
    scale_color_manual(values = c("grey10", "black")) +
    labs(title = "",
         x = "",
         y = "Maximal respiration \n(pmol/min per worm)") +
    theme_classic() +
    scale_x_discrete(labels = c("N2", "*swip-10*")) +
    theme(axis.text.x = ggtext::element_markdown(),
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 11),
          legend.position = "none")
#maximal
#dev.off()

# SPARE CAPACITY
#png("swip-10/figures/spare.png", width = 3, height = 4, units = "in", res = 300)
spare <- respiratory.measures %>% 
    filter(label != "swip20") %>% 
    ggplot(mapping  = aes(x= group, y  = spare))  +
    geom_boxplot(aes(fill  = group, color = group), width = 0.5)  +
    geom_jitter(aes(color = group), shape  = 19, width = 0.1) +
    scale_fill_manual(values = c("black", "#8B0000")) +
    scale_color_manual(values = c("grey10", "black")) +
    labs(title = "",
         x = "",
         y = "Spare capacity \n(pmol/min per worm)"
    ) +
    theme_classic() +
    scale_x_discrete(labels = c("N2", "*swip-10*")) +
    theme(axis.text.x = ggtext::element_markdown(),
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 11),
          legend.position = "none")
#spare
#dev.off()

### Table with raw data
ocr_raw %>% 
    filter(label != "BG") %>% 
    write_tsv(.,"raw_files/for_dryad/seahorse_raw.txt")

respiratory.measures %>% 
    write_tsv(.,"raw_files/for_dryad/seahorse_respiratory_measures.txt")

png("swip-10/figures/seahorse_all_arranged.png", width = 5.5, height = 5.5, res = 300, units = "in")

bottom <- plot_grid(basal, maximal, spare, nrow = 1, labels = c("B", "C", "D"), label_size = 12)

plot_grid(profile, bottom, ncol = 1, labels = c("A", ""), label_size = 12)

dev.off()
