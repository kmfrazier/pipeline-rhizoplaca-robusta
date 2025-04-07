library(ggplot2)
library(readr)

data <- read_csv("fig1.csv")

# Open a new window for the plot
x11()

g <- ggplot(data, aes(x = gene, y = log(`E-value`, base = exp(1))), show.legend = TRUE) + # nolint
  geom_point(color = "blue", size = 1) +
  geom_hline(yintercept = -5, color = "blue", linetype = "dotted", show.legend = TRUE, size = 1) + # nolint
  labs(title = "tBLASTn Hits of Trebouxia aggregata Photosystem I Genes", x = "Gene Name", y = "Natural Log of E-value") + # nolint
  annotate("text", x = 1.5, y = -4.3, label = "Significant Match Threshhold", color = "blue") + # nolint
  theme(legend.position = "right")

plot(g)
