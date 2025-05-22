library(dplyr)
library(highcharter)
library(htmlwidgets)
library(webshot)


df <- read.csv("D:/BaiduNetdiskDownload/cvd_count.csv")

df_filtered <- df[!grepl("\\.", df$id), ]
df_filtered$auc_value <- as.numeric(sub(" .*", "", df_filtered$AUC))
bins <- c(0, 0.6, 0.7, 0.8, 0.9, 1.0)
labels <- c("<0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "≥0.9")
df_filtered$range <- cut(df_filtered$auc_value, 
                         breaks = bins, 
                         labels = labels, 
                         right = FALSE, 
                         include.lowest = TRUE)
range_counts <- table(df_filtered$range)


all_labels <- data.frame(label = labels)
range_counts_df <- data.frame(label = names(range_counts), count = as.vector(range_counts))
df_percent <- all_labels %>%
  left_join(range_counts_df, by = "label") %>%
  mutate(count = ifelse(is.na(count), 0, count))


total_count <- sum(df_percent$count)


if (total_count == 0) {
  actual_percentages <- rep(0, nrow(df_percent))
} else {
  actual_percentages <- (df_percent$count / total_count) * 100
}


df_percent$actual_percentage <- actual_percentages
df_percent$whole_percentage <- floor(actual_percentages)
df_percent$fractional_part <- actual_percentages - floor(actual_percentages)


sum_whole <- sum(df_percent$whole_percentage)
difference <- 100 - sum_whole


if (difference > 0) {
  sorted_indices <- order(df_percent$fractional_part, decreasing = TRUE)
  indices_to_increase <- sorted_indices[1:difference]
  df_percent$whole_percentage[indices_to_increase] <- df_percent$whole_percentage[indices_to_increase] + 1
}

percentages <- df_percent$whole_percentage
data_pair <- data.frame(label = labels, percentage = percentages)


range_meanings <- lapply(labels, function(label) {
  meanings <- df_filtered$meaning[df_filtered$range == label]
  paste(meanings, collapse = "<br>")
})
names(range_meanings) <- labels


colors <- c("#CE93D8", "#FFCC80", "#82B1FF", "#80DEEA", "#F48FB1")


data <- list()
for (i in 1:nrow(data_pair)) {
  data[[i]] <- list(
    name = data_pair$label[i],
    y = data_pair$percentage[i],
    meanings = range_meanings[[data_pair$label[i]]],
    color = colors[i]
  )
}


hc <- highchart() %>%
  hc_chart(type = "pie", width = 1200, height = 600) %>%  
  hc_plotOptions(
    pie = list(
      allowPointSelect = TRUE,
      showInLegend = TRUE,
      legendType = "point"
    )
  ) %>%
  hc_add_series(
    name = "AUC Ranges",
    data = data,
    innerSize = "0%",
    size = "50%",
    showInLegend = TRUE,
    dataLabels = list(
      enabled = TRUE,
      useHTML = TRUE,        
      crop = FALSE,         
      overflow = "allow",
      format = "<div style='white-space: nowrap; text-align: left;'>{point.meanings}</div>",
      distance = 28,
      style = list(fontSize = "10px", fontFamily = "Arial", fontWeight = "normal"),  
      backgroundColor = "#f8f8f8",
      borderWidth = 0.1,
      borderRadius = 1,
      padding = 1
    )
  ) %>%
  hc_add_series(
    name = "AUC Percentages",
    data = data,
    innerSize = "0%",
    size = "50%",
    showInLegend = FALSE,
    dataLabels = list(
      enabled = TRUE,
      format = "{point.percentage:.0f}%",
      distance = -20,
      style = list(fontSize = "12px", fontFamily = "Arial", fontWeight = "normal")  
    )
  ) %>%
  hc_legend(
    enabled = TRUE,
    align = "center",
    verticalAlign = "bottom",
    layout = "horizontal",
    # y = -180,  
    labelFormatter = htmlwidgets::JS("function(){ return this.name; }"),
    itemMarginTop = 2,
    itemMarginBottom = 2,
    itemStyle = list(fontSize = "10px", fontFamily = "Arial", fontWeight = "normal"),
    symbolHeight = 12,
    symbolWidth = 12,
    symbolRadius = 0,
    itemDistance = 6   
  ) %>%
  hc_title(text = "") %>%
  hc_exporting(
    enabled = TRUE,
    filename = "auc_pie",
    scale = 9,  
    buttons = list(
      contextButton = list(
        menuItems = list("downloadSVG", "downloadPNG", "downloadJPEG")
      )
    )
  )

hc

htmlwidgets::saveWidget(hc, "pie.html", selfcontained = TRUE)
browseURL("pie.html")



