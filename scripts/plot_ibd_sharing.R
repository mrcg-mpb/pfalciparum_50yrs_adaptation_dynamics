

library(tidyverse)

col.code <- c("#000000", "#FFCC00") 

ibd <- read_delim("from_alfred/gambia1966_72_ibdSharing_majAllele.tsv")

ibd <- ibd %>% 
   # Compute chromosome size
   group_by(chr) %>% 
   summarise(chr_len = max(pos_bp)) %>% 
   
   # Calculate cumulative position of each chromosome
   mutate(tot = cumsum(chr_len) - chr_len) %>%
   select(-chr_len) %>%
   
   # Add this info to the initial dataset
   left_join(ibd, ., by=c("chr" = "chr")) %>%
   
   # Add a cumulative position of each SNP
   arrange(chr, pos_bp) %>%
   mutate( BPcum = pos_bp + tot, 
           chr = as.integer(factor(chr))) %>% 
   select(-tot)

axis <- ibd %>%
   group_by(chr) %>%
   summarize(center = mean(BPcum))

ylimits <- c(floor(min(ibd$log10_pvalue)), floor(max(ibd$log10_pvalue)) + 1)


ibd %>% 
   ggplot(aes(x=BPcum, y=log10_pvalue, color = as_factor(chr))) +
   # Show all points
   geom_point(size = 2) + # alpha = 0.75, 

   # custom X axis:
   scale_x_continuous(label = axis$chr, breaks = axis$center) +
   scale_y_continuous(expand = c(0,0), limits = ylimits) +
   scale_color_manual(values = rep(col.code, unique(length(axis$chr)))) +
   # scale_size_continuous(range = c(0.5,3)) +
   labs(x = "", y = expression(-~log10("P-values"))) + 
   
   # Custom the theme:
   theme_classic() +
   theme(legend.position = "none",
         axis.text = element_text(size = 9, color = "#000000"),
         axis.title = element_text(face = "bold"),
         axis.line = element_line(color = "#000000"))

ggsave("from_alfred/gambia1966_IBDSharing_majAllele.pdf", 
       units = "mm", width = 190, height = 150, dpi = 600)
