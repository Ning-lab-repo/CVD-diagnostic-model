
library(igraph)
library(ggraph)
library(ggplot2)



data <- read.csv("D:/BaiduNetdiskDownload/cvd_count.csv", stringsAsFactors = FALSE)


second_layer <- data[!grepl("\\.", data$id), ]
third_layer <- data[grepl("\\.", data$id), ]


root <- "CVDs"


edges_root_to_second <- data.frame(from = rep(root, nrow(second_layer)), 
                                   to = second_layer$id, 
                                   stringsAsFactors = FALSE)


third_layer$prefix <- sub("\\..*", "", third_layer$id)


edges_second_to_third <- data.frame(from = third_layer$prefix, 
                                    to = third_layer$id, 
                                    stringsAsFactors = FALSE)


edges_second_to_third <- edges_second_to_third[edges_second_to_third$from %in% second_layer$id, ]


edges <- rbind(edges_root_to_second, edges_second_to_third)


root_label <- ""
second_labels <- paste(second_layer$Disease, " (", second_layer$Count, ")", sep = "")
third_labels <- paste(third_layer$Disease, " (", third_layer$Count, ")", sep = "")

node_data <- data.frame(id = c(root, second_layer$id, third_layer$id),
                        label = c(root_label, second_labels, third_labels),
                        stringsAsFactors = FALSE)


g <- graph_from_data_frame(edges, directed = TRUE, vertices = node_data)


layout <- layout_as_tree(g, root = root)

root_index <- which(node_data$id == root)
second_indices <- which(node_data$id %in% second_layer$id)
third_indices <- which(node_data$id %in% third_layer$id)

compress_factors <- c(
  root_to_second = 0.1,   
  second_to_third = 1     
)


original_root_height <- layout[root_index, 2]
layout[root_index, 2] <- original_root_height * compress_factors["root_to_second"]


layout[second_indices, 2] <- layout[second_indices, 2] * compress_factors["root_to_second"]


for (node in third_layer$id) {
  parent <- edges$from[edges$to == node]
  parent_y <- layout[which(node_data$id == parent), 2]
  current_idx <- which(node_data$id == node)
  offset <- (layout[current_idx, 2] - parent_y) * compress_factors["second_to_third"]
  layout[current_idx, 2] <- parent_y + offset
}


new_layout <- layout
#new_layout[,1] <- -layout[,2]  # x = -y
#new_layout[,2] <- -layout[,1]   # y = x


horizontal_scale <- 0.6  
new_layout[,1] <- new_layout[,1] * horizontal_scale


new_layout[,2] <- new_layout[,2] * 0.9


p <- ggraph(g, layout = new_layout) +
  geom_edge_diagonal(
    color = "grey",
    alpha = 0.6,
    width = 0.8  
  ) +
  geom_node_point(
    color = "black", 
    fill = "violetred3", 
    shape = 21, 
    size = 4, 
    stroke = 0.8
  ) +
  geom_node_text(
    aes(label = label),
    hjust = 0,
    nudge_x = 0, 
    nudge_y = 0.002,    
    size = 4,
    angle = 0,
    check_overlap = TRUE
  ) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = unit(c(1,14, 1, 1), "cm")) + 
  scale_x_reverse() +
  scale_y_reverse() +  
  coord_flip()  

print(p)


ggsave("D:/BaiduNetdiskDownload/tree.pdf", p, width = 20, height = 12, device = "pdf")
