


plot_betaScore <- function(data,
                          format = "pdf",
                          ...
                        )
{

   data <- data %>%
      # Compute chromosome size
      group_by(Chr) %>%
      summarise(chr_len = max(Start)) %>%

      # Calculate cumulative position of each chromosome
      mutate(tot = cumsum(chr_len) - chr_len) %>%
      select(-chr_len) %>%

      # Add this info to the initial dataset
      left_join(data, ., by=c("Chr" = "Chr")) %>%

      # Add a cumulative position of each SNP
      arrange(Chr, Start) %>%
      mutate( BPcum = Start + tot)

   axis_set <- data %>%
      group_by(Chr) %>%
      summarize(center = mean(BPcum))

   ylimits <- c(floor(min(data$betaScore)), floor(max(data$betaScore)) + 1)

   col.code <- c("gray40", "#0099FF")

   data %>%
      ggplot(aes(x=BPcum, y=betaScore, color = as_factor(Chr))) +
      # Show all points
      geom_point(size = 1.2) + # alpha = 0.75,
      geom_hline(yintercept = 0, color = "grey40", linetype = "dashed") +
      # custom X axis:
      scale_x_continuous(label = axis_set$Chr, breaks = axis_set$center) +
      # remove space between plot area and x axis
      scale_y_continuous(expand = c(0,0), limits = ylimits) +
      scale_color_manual(values = rep(col.code, unique(length(axis_set$Chr)))) +
      scale_size_continuous(range = c(0.5,3)) +
      labs(x = "", y = expression("Beta Scores" ~ beta ~ "1*")) +
      # Custom the theme:
      theme_classic() +
      theme(legend.position = "none",
            axis.text = element_text(size = 8),
            axis.title = element_text(face = "bold"))
}
