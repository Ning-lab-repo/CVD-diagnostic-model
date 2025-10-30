
library(stringr)
library(pheatmap)
library(dplyr)



df_sample <- read.csv("D:/BaiduNetdiskDownload/ukbnmr_325_data.csv", header = TRUE, stringsAsFactors = FALSE)


nmr_group <- read.csv("D:/BaiduNetdiskDownload/heatmap_group_ukbnmr.csv", header = TRUE, stringsAsFactors = FALSE)

colnames(nmr_group) <- c("Metabolite", "Group")


desired_groups <- c("Amino acids", "Apolipoproteins", "Cholesterol", "Cholesteryl esters", "Free cholesterol", "Fatty acids", "Fluid balance", 
                    "Glycolysis related metabolites", "Inflammation", "Ketone bodies", "Lipoprotein particle concentrations", "Lipoprotein particle sizes",
                    "Lipoprotein subclasses", "Relative lipoprotein lipid concentrations", "Other lipids", "Phospholipids", "Total lipids", "Triglycerides", 
                    "Relative lipid concentrations", "Relative cholesterol concentrations", "Relative lipoprotein cholesterol concentrations")

nmr_group$Group <- factor(nmr_group$Group, levels = desired_order)

nmr_group_sorted <- nmr_group %>%
  arrange(Group, Metabolite) %>%
  filter(Metabolite %in% colnames(df_sample))

ordered_metabolites <- nmr_group_sorted$Metabolite


#df_sample$illne_reported <- as.character(df_sample$illne_reported)
df_sample$diagnosis_cvd_id <- as.character(df_sample$diagnosis_cvd_id)


excluded_codes <- c('I79','I80','I81','I82','I83','I84',
                    'I85','I86','I87','I88','I89','I95','I97','I98','I99')


df_sample$diagnosis <- apply(df_sample, 1, function(row) {

  diagnoses <- c(strsplit(as.character(row['diagnosis_cvd_id']), '\\|')[[1]])
  

  has_valid_cvd <- any(sapply(diagnoses, function(diagnosis) {
  
    (grepl("^I", diagnosis) && !any(sapply(excluded_codes, function(code) grepl(paste("^", code, sep=""), diagnosis)))) || grepl("^G45", diagnosis)
  }))
  
  if (has_valid_cvd) {
    return(1)
  } else {
    return(0)
  }
})

cvd_count <- sum(df_sample$diagnosis == 1, na.rm = TRUE)
print(paste("Diagnosis:", cvd_count))


disease_codes <- c('I10','I25','I20','I48','I21','I50','I47','I73','I26','G45','I51','I63','I49','I44','I34','I35','I45','I67','I77','I12','I64','I42','I69','I60','I70','I22',
                   'I65','I74','I24','I31','I08','I71','I61','I78','I46','I27','I05')


disease_groups <- lapply(disease_codes, function(code) {
  apply(df_sample, 1, function(row) {
    check_code(row['diagnosis_cvd_id'], code)
  })
})
names(disease_groups) <- disease_codes


set.seed(123)
sample_group <- function(indices, max_samples = 100) {
  if(length(indices) > max_samples) sample(indices, max_samples) else indices
}


sampled_indices <- list()


non_cvd_indices <- which(df_sample$diagnosis == 0)
if(length(non_cvd_indices) > 0) {
  sampled_indices[["non_CVD"]] <- sample_group(non_cvd_indices)
}

for(disease in names(disease_groups)) {
  indices <- which(disease_groups[[disease]])
  if(length(indices) > 0) {
    sampled_indices[[disease]] <- sample_group(indices)
  }
}


all_samples <- unlist(sampled_indices)
original_ids <- rownames(df_sample)[all_samples]


make_unique_ids <- function(ids) {
  ave(ids, ids, FUN = function(x) {
    if(length(x) > 1) paste0(x, "_", seq_along(x)) else x
  })
}

unique_ids <- make_unique_ids(original_ids)


metabolite_data <- df_sample[all_samples, ordered_metabolites]
zscored <- scale(metabolite_data)
#zscored[zscored > 2] <- 2
#zscored[zscored < -2] <- -2
rownames(zscored) <- unique_ids





library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(ggsci)
library(randomcoloR)

col_fun <- colorRamp2(
  breaks = c(-2, 0, 2), 
  colors = c("blue", "white", "red")
)


group_colors <- c(
  "Cholesterol" = "#FF77FF",
  "Cholesteryl esters" = "#FF3A7C",
  "Free cholesterol" = "#80DEEA",
  "Phospholipids" = "#26A69A",
  "Lipoprotein particle concentrations" = "#82B1FF",
  "Total lipids" = "#FFA000",
  "Other lipids" = "#FFD54F",
  "Apolipoproteins" = "#8E24AA",
  "Lipoprotein subclasses" = "#7CB342",
  "Amino acids" = "#FF6B6B",
  "Fatty acids" = "#4ECDC4",
  "Fluid balance" = "#45B7D1",
  "Glycolysis related metabolites" = "#96CEB4",
  "Inflammation" = "#FFEEAD",
  "Ketone bodies" = "#D4A5A5",
  "Lipoprotein particle sizes" = "#9B59B6",
  "Relative lipoprotein lipid concentrations" = "#3498DB",
  "Triglycerides" = "#E74C3C",
  "Relative lipid concentrations" ="#F48FB1",
  "Relative cholesterol concentrations" = "#DCE775",
  "Relative lipoprotein cholesterol concentrations" = "#B39DDB"
  
)

group_labels <- rep(
  names(sampled_indices),
  times = sapply(sampled_indices, length)
)


levels_new <- c("non_CVD", setdiff(unique(group_labels), "non_CVD"))
group_labels <- factor(group_labels, levels = levels_new)
ord <- order(group_labels)  
mat_reordered <- t(zscored)[, ord]
group_labels_reordered <- group_labels[ord]


num_diseases <- length(unique(group_labels_reordered))
print(paste("Number of unique diseases:", num_diseases))

#base_colors <- colorRampPalette(pal_npg("nrc")(10))(num_diseases)
#base_colors <- distinctColorPalette(num_diseases)
#disease_colors <- setNames(colorRampPalette(base_colors)(num_diseases), unique(group_labels_reordered))
# Disease colors
num_diseases <- length(unique(group_labels_reordered))
base_colors <- distinctColorPalette(num_diseases)
disease_colors <- setNames(colorRampPalette(base_colors)(num_diseases), unique(group_labels_reordered))
# Override non_CVD color specifically
disease_colors["non_CVD"] <- "#0040FF"

row_anno <- rowAnnotation(
  `Biomarker group` = nmr_group_sorted$Group,  
  col = list(`Biomarker group` = group_colors),  
  show_annotation_name = TRUE,  
  annotation_name_side = "bottom",  
  annotation_name_gp = gpar(fontsize = 14, family = "Arial"),  
  annotation_legend_param = list(
    title = "Biomarker group", 
    title_gp = gpar(fontsize = 16, fontface = "bold", family = "Arial"),
    labels_gp = gpar(fontsize = 14, family = "Arial"),
    title_position = "topleft",
    direction = "vertical"
  )
)


col_anno <- HeatmapAnnotation(
  df = data.frame(Disease = group_labels_reordered),
  col = list(Disease = disease_colors),
  show_annotation_name = TRUE,  
  annotation_name_side = "right", 
  annotation_name_gp = gpar(fontsize = 14, family = "Arial"),  
  annotation_legend_param = list(
    title = "Disease",
    title_gp = gpar(fontsize = 16, fontface = "bold", family = "Arial"),
    labels_gp = gpar(fontsize = 14, family = "Arial"),
    title_position = "topleft",
    direction = "vertical"
  )
)


pdf("Heatmap_Sampling_100.pdf", width = 48, height = 78) 

ht <- Heatmap(
  matrix = mat_reordered,
  name = "Z-score",
  col = col_fun,
  
  

  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,      
  row_names_side = "right",      
  column_title = NULL,           
  

  row_split = nmr_group_sorted$Group,
  row_title = NULL,            
  column_split = group_labels_reordered,
  
  
  top_annotation = col_anno,      
  left_annotation = row_anno,     
  
  row_names_gp = gpar(fontsize = 16, family = "Arial"),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 16, fontface = "bold", family = "Arial"),
    labels_gp = gpar(fontsize = 14, family = "Arial"),
    legend_height = unit(4, "cm"),
    title_position = "topleft"
  )
)


draw(ht, 
     padding = unit(c(5, 5, 5, 5), "cm"),
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     merge_legends = TRUE,           
     auto_adjust = FALSE,            
     legend_grouping = "adjusted"
)

dev.off()
