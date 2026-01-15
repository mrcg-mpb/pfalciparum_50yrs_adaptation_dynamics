

## Load packages
library(tidyverse)
library(stringr)
library(ggrepel)
library(glue)
library(rlang)


# load New Gene Name file
Gene_name <- "data/PF_GeneName.tsv"
Gene_name <- read_tsv(Gene_name)

## Let's combined all results in on single file
## Single country population

## I can also use map_dfr to combine all dataframes into one

# Data <- readr::read_tsv("keneba66_72_fuli.xlsx")

Data <- Data  %>% 
   
   # Compute chromosome size
   group_by(Chrom) %>% 
   summarise(chr_len = max(Start)) %>% 
   
   # Calculate cumulative position of each chromosome
   mutate(tot = cumsum(chr_len) - chr_len) %>%
   select(-chr_len) %>%
   
   # Add this info to the initial dataset
   left_join(Data, ., by=c("Chrom" = "Chrom"))  %>% 
   
   # Add a cumulative position of each SNP
   arrange(Chrom, Start) %>%
   mutate( BPcum = Start + tot) %>% 
   
   # Add gene name
   dplyr::select(-c(tot, Chrom))

# Then we need to prepare the X axis. Indeed we do not want to display the cumulative position of SNP in bp, 
# but just show the chromosome name instead.

axis_set <- Data %>%
   group_by(Chr) %>%
   summarize(center = mean(BPcum))

ylimits <- c(floor(min(Data$TajimasD)), floor(max(Data$TajimasD)) + 1)

# Ready to make the manhattan plot using ggplot2:
tajima.manhattan <- function(data, xcolumn, ycolumn, color, axis_set, ylimits)
{
   data %>% 
      ggplot(aes(x = {{xcolumn}}, y = {{ycolumn}}, color = as_factor({{color}}))) + 
      
      # Show all points
      geom_point(alpha = 1.2) +
      geom_hline(yintercept = 0, color = "black", linetype = "dashed") + 
      
      # custom X axis:
      scale_x_continuous(label = axis_set$Chr, breaks = axis_set$center) +
      scale_y_continuous(expand = c(0,0), limits = {{ylimits}}) + # remove space between plot area and x axis
      scale_color_manual(values = rep(c("#000000", "#FFCC00") , unique(length(axis_set$Chr)))) +
      scale_size_continuous(range = c(0.5,3)) +
      labs(x = "Chromosomes", y = "Tajima's D") + 
      
      # Custom the theme:
      theme_classic() +
      theme(legend.position = "none",
            axis.text = element_text(face = "bold"),
            # panel.grid.major.x = element_blank(),
            # panel.grid.minor.x = element_blank(),
            axis.title = element_text(face = "bold.italic"))
}

tajima.manhattan(Data, BPcum, TajimasD, Chr, axis_set, ylimits)

ggsave("from_alfred/gambia1966_72_tajima.pdf", units = "mm",
       width = 200, height = 100, dpi = 600)

# Explore top genes
threshold <- abs(2*mean(Data$TajimasD) + sd(Data$TajimasD))
Data %>% 
   dplyr::filter(TajimasD >= threshold) %>% 
   as_tibble() %>% 
   print(n=50)

threshold <- quantile(Data$TajimasD, 0.95)
Data %>% 
   dplyr::filter(TajimasD >= threshold) %>% 
   as_tibble() %>% 
   print(n=50)

#==========================================
#==========================================
# VERSION 1

Data <- read_tsv("gambia1966_tajima_v1.tsv") %>% 
   dplyr::mutate(Gene = gsub("SNPEFF_TRANSCRIPT_ID=", "", Gene),
                 Gene = sub("\\.1", "", Gene)) %>% 
   dplyr::inner_join(., Gene_name, by = c("Gene" = "Gene_ID"))

Data <- Data  %>% 
   
   # Compute chromosome size
   group_by(Chromosome) %>% 
   summarise(chr_len = max(Start)) %>% 
   
   # Calculate cumulative position of each chromosome
   mutate(tot = cumsum(chr_len) - chr_len) %>%
   select(-chr_len) %>%
   
   # Add this info to the initial dataset
   left_join(Data, ., by=c("Chromosome" = "Chromosome"))  %>% 
   
   # Add a cumulative position of each SNP
   arrange(Chromosome, Start) %>%
   mutate( BPcum = Start + tot) %>% 
   rename(TajimasD = tajimasd) %>% 
   
   # Remove unneeded columns
   select(-c(H, S, theta, khat, tot, Chromosome, NewGeneName, `#Iso_used`, `Total#Iso`, `Iso%`))

# Then we need to prepare the X axis. Indeed we do not want to display the cumulative position of SNP in bp, 
# but just show the chromosome name instead.

axis_set <- Data %>%
   group_by(Chr) %>%
   summarize(center = mean(BPcum))

ylim <- floor(max(Data$TajimasD)) + 1
ylimits <- c(floor(min(Data$TajimasD)), ylim)

tajima.manhattan(Data, BPcum, TajimasD, Chr, axis_set, ylimits)

ggsave("gambia1966_tajima_v1.pdf", units = "mm",
       width = 190, height = 150, dpi = 600)







