library(ggplot2)
library(dplyr)
library(readr)


data <- read_csv("fig3.csv")

data_count <- data %>%
  group_by(gene_name, expect) %>%
  summarize(count = n(), .groups = "drop")

x11()

g <- ggplot(data_count, aes(x = gene_name, y = log(as.numeric(expect), base = exp(1)), size = count, color = gene_name)) + # nolint
  geom_point() +
  geom_hline(yintercept = -5, color = "blue", linetype = "dotted", size = 1) +
  labs(
    title = "BLASTx Hits of Trebouxia aggregata Photosynthesis Related Genes",
    x = "Gene Name",
    y = "Natural Log of Expect Value",
    color = "Gene Name",
    size = "Number of Points"
  ) +
  annotate("text", x = 4, y = -4.3, label = "Significant Match Threshold", color = "blue") + # nolint
  theme(legend.position = "right") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

plot(g)
