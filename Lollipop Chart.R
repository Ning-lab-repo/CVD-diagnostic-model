
library(ggplot2)


data <- read.csv("D:/BaiduNetdiskDownload/cvd_count.csv", check.names = FALSE)


data$AUC_value <- as.numeric(sub(" \\(.*\\)", "", data$AUC))

filtered_data <- data[!grepl("\\.", data$id), ]


filtered_data <- filtered_data[order(filtered_data$AUC_value), ]


filtered_data$meaning <- factor(filtered_data$meaning, levels = filtered_data$meaning)


p <- ggplot(filtered_data, aes(x = AUC_value, y = meaning)) +
  geom_segment(aes(x = 0.5, xend = AUC_value, y = meaning, yend = meaning, color = AUC_value), size = 1) +
  geom_point(aes(color = AUC_value), size = 3) +
  scale_color_gradient(low = "#BFDFF5", high = "#608CB9") +
  geom_text(aes(label = sprintf("%.4f", AUC_value)), hjust = -0.2) +
  scale_x_continuous(limits = c(0.5, 1)) +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "grey95", size = 0.2),
    panel.grid.minor = element_line(color = "grey95", size = 0.1),
    legend.position = "right",
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.text.y = element_text(margin = margin(r = -10)) 
  ) +
  labs(x = "", y = "", color = "AUC")

print(p)
ggsave("D:/BaiduNetdiskDownload/Lollipop Chart.pdf", p, width = 10, height = 8, device = "pdf")


