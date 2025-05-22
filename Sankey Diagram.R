library(networkD3)
library(dplyr)
library(readr)
library(tidyr)
library(webshot)

create_sorted_sankey_diagram <- function(csv_file) {

  df <- read_csv(csv_file)
  

  
  
  left_cols <- c('IDL_CE','LDL_C_pct','LA','IDL_C','LA_pct','L_LDL_C_pct','XS_VLDL_FC_pct_C','XS_VLDL_CE_pct','IDL_CE_pct','XS_VLDL_C_pct','XS_VLDL_CE_pct_C','Creatinine','M_LDL_C_pct')
  
  middle_ids <- c(
    'Multiple valve diseases',
    'Essential (primary) hypertension',
    'Hypertensive renal disease',
    'Angina pectoris',
    'Acute myocardial infarction',
    'Subsequent myocardial infarction',
    'Other acute ischaemic heart diseases',
    'Chronic ischaemic heart disease',
    'Nonrheumatic aortic valve disorders',
    'Cardiomyopathy',
    'Cardiac arrest',
    'Atrial fibrillation and flutter',
    'Heart failure',
    'Complications and ill-defined descriptions of heart disease',
    'Cerebral infarction',
    'Stroke, not specified as haemorrhage or infarction',
    'Occlusion and stenosis of precerebral arteries, not resulting in cerebral infarction',
    'Atherosclerosis',
    'Aortic aneurysm and dissection',
    'Other peripheral vascular diseases',
    'Other disorders of arteries and arterioles'
  )
  
  #right_cols <- c('T-Bil', 'RDW', 'Chol', 'BMI', 'Age', 'DBil', 'WC', 'Sex', 'Lymph%', 'Neutro%', 'MCV', 'TG', 'Neutro-C', 'WBC-C')
  right_cols <- c('CYS','Age','GGT','HBA1C','RDW','WC','BUN','CHOL','TES','LDLD','MO')

  left_cols <- left_cols[left_cols %in% colnames(df)]
  right_cols <- right_cols[right_cols %in% colnames(df)]
  middle_ids <- middle_ids[middle_ids %in% df$id]
  

  links <- data.frame(source = character(), target = character(), value = numeric(), stringsAsFactors = FALSE)
  

  for (left_col in left_cols) {
    for (mid_id in middle_ids) {
      weight <- df %>% filter(id == mid_id) %>% pull(!!sym(left_col))
      if (!is.na(weight) && weight != 0) {
        links <- rbind(links, data.frame(source = left_col, target = mid_id, value = abs(weight)))
      }
    }
  }
  

  for (mid_id in middle_ids) {
    for (right_col in right_cols) {
      weight <- df %>% filter(id == mid_id) %>% pull(!!sym(right_col))
      if (!is.na(weight) && weight != 0) {
        links <- rbind(links, data.frame(source = mid_id, target = right_col, value = abs(weight)))
      }
    }
  }
  

  left_flows <- links %>%
    filter(source %in% left_cols) %>%
    group_by(source) %>%
    summarise(total_flow = sum(value)) %>%
    arrange(desc(total_flow))

  middle_flows <- links %>% 
    filter(source %in% left_cols, target %in% middle_ids) %>%
    group_by(target) %>%
    summarise(total_flow = sum(value)) %>%
    arrange(desc(total_flow))
  print(middle_flows)
  
  right_flows <- links %>%
    filter(target %in% right_cols) %>%
    group_by(target) %>%
    summarise(total_flow = sum(value)) %>%
    arrange(desc(total_flow))
  

  n_left <- nrow(left_flows)
  n_middle <- nrow(middle_flows)
  n_right <- nrow(right_flows)
  
  
  nodes <- data.frame(name = c(left_flows$source, middle_flows$target, right_flows$target))
  

  nodes$group <- nodes$name
  

  light_colors_extended <- c(
    "#FFAB91", # Light Pink
    "#FFF9C4", # Light Yellow
    "#F48FB1", # Light Blue
    "#FFB3A7", # Light Green
    "#E1BEE7", # Light Orange
    "#B3E5FC",  # Light Grey Blue
 #   "#B2EBF2", # Light Cyan Blue
    "#F8BBD0", # Light Salmon
    "#B39DDB", # Light Lemon
    "#EDE7F6", # Light Lavender
    "#E0E0E0", # Light Silver
    "#F0F4C3", # Light Lime
    "#DCE775", # Lime 2
    "#C8E6C9", # Light Mint
    "#B2DFDB", # Light Teal
    "#A5D6A7", # Light Olive
    "#90CAF9", # Light Sky Blue
    "#81D4FA", # Light Azure
    "#80CBC4", # Teal 2
    "#80DEEA", # Cyan 2
    "#82B1FF", # Light Indigo
    "#9FA8DA", # Light Steel Blue
    "#CE93D8", # Light Violet
    "#E1BEE7", # Light Mauve
    "#FFAB91", # Light Coral
    #  "#FFCC80", # Light Amber
    #  "#FFECB3", # Light Gold
    "#FFE0B2" # Light Peach
    #  "#FFD54F", # Light Yellow 2
    #  "#FFF59D" # light Yellow 3
  )
  
  color_scale <- paste0("d3.scaleOrdinal().range(['", paste(light_colors_extended, collapse = "','"), "'])")
  #color_scale <- 'd3.scaleOrdinal(d3.schemePastel1)'

  links$source <- match(links$source, nodes$name) - 1
  links$target <- match(links$target, nodes$name) - 1
  links$group <- nodes$group[links$source + 1]
  

  sankey <- sankeyNetwork(
    Links = links,
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "name",
    NodeGroup = "group",
    LinkGroup = "group",
    colourScale = JS(color_scale),
    fontFamily = "Arial",
    nodeWidth = 20,
    fontSize = 18,
    width = 1500,   
    height = 900,   
    nodePadding = 16,
    #height = target_height,
    iterations = 0
  )
  
  return(sankey)
}


graph <- create_sorted_sankey_diagram("D:/BaiduNetdiskDownload/sankey_data.csv")
graph

saveNetwork(graph, "sankey.html")

webshot::webshot("sankey.html", "sankey.png", vwidth = 1500, vheight = 600)
