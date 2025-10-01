# Libraries
#==========================================================================================================
library(ggplot2)
library(ggstatsplot)
library(dplyr)
library(ggdist)
library(reshape2)

#Functions
#==========================================================================================================

# Get pathway statistics
load_rdata_from_subdirectories <- function(main_dir) {

  sub_dirs <- list.dirs(main_dir, full.names = TRUE, recursive = FALSE)

  for (sub_dir in sub_dirs) {
    rdata_files <- list.files(sub_dir, pattern = "\\.RData$", full.names = TRUE)

    for (file in rdata_files) {
      load(file, envir = .GlobalEnv)  # Load the .RData file into the environment
      cat("Loaded:", file, "\n")
    }
  }
}

# Plot statistics
plot_statistics <- function(
    files,
    group_labels,
    levels_order = NULL,
    output_file = NULL,
    colors = NULL,
    width = 7,
    height = 5
) {
  # Load files
  lapply(files, load, envir = .GlobalEnv)

  # Get score metric variables
  median_values <- lapply(files, function(f) {
    obj_name <- basename(f) %>%
      sub("^statistics_", "", .) %>%
      sub("(_GSE[0-9]+)?\\.RData$", "", .)
    get(obj_name)$Median
  })

  data <- data.frame(
    Median = unlist(median_values),
    Group = rep(group_labels, times = sapply(median_values, length))
  )

  if (!is.null(levels_order)) {
    data$Group <- factor(data$Group, levels = levels_order)
  }

  # Plot
  p <- ggbetweenstats(
    data = data,
    x = Group,
    y = Median,
    type = "nonparametric",
    plot.type = "box",
    pairwise.comparisons = FALSE,
    pairwise.display = "significant",
    centrality.plotting = FALSE,
    bf.message = FALSE,
    point.args = list(
      position = ggplot2::position_jitterdodge(dodge.width = 0.6, jitter.width = 0.1),
      alpha = 0.8,
      size = 1,
      stroke = 0,
      na.rm = TRUE
    ),
    boxplot.args = list(
      width = 0.1,
      alpha = 0.2,
      na.rm = TRUE
    )
  ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 18),
      axis.text.y = ggplot2::element_text(size = 18),
      plot.subtitle = ggplot2::element_text(size = 12)
    )

  if (!is.null(colors)) {
    p <- p + ggplot2::scale_color_manual(values = colors)
  }

  # Save or return
  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, plot = p, width = width, height = height)
  } else {
    return(p)
  }
}


# Figure 2A
#==========================================================================================================
# For GSE48080
plot_statistics(
  files = c(
    "data/KEGG/statistics_score_metrics_survivors_1_7_GSE48080.RData",
    "data/KEGG/statistics_score_metrics_non-survivors_1_7_GSE48080.RData"
  ),
  group_labels = c("S 1-7", "NS 1-7"),
  levels_order = c("S 1-7", "NS 1-7"),
  output_file = "GSE48080_KEGG_stat.svg",
  colors = c("snow4", "violetred")
)

plot_statistics(
  files = c(
    "data/GO/statistics_score_metrics_survivors_1_7_GSE48080.RData",
    "data/GO/statistics_score_metrics_non-survivors_1_7_GSE48080.RData"
  ),
  group_labels = c("S 1-7", "NS 1-7"),
  levels_order = c("S 1-7", "NS 1-7"),
  output_file = "GSE48080_GO_stat.svg",
  colors = c("snow4", "violetred")
)

# For GSE54514
plot_statistics(
  files = c(
    "data/KEGG/statistics_score_metrics_control_1_5_GSE54514.RData",
    "data/KEGG/statistics_score_metrics_survivors_1_5_GSE54514.RData",
    "data/KEGG/statistics_score_metrics_non-survivors_1_5_GSE54514.RData"
  ),
  group_labels = c("HC 1-5", "S 1-5", "NS 1-5"),
  levels_order = c("HC 1-5", "S 1-5", "NS 1-5"),
  output_file = "GSE54514_KEGG_stat.svg",
  colors = c("yellowgreen", "snow4", "violetred")
)

plot_statistics(
  files = c(
    "data/GO/statistics_score_metrics_control_1_5_GSE54514.RData",
    "data/GO/statistics_score_metrics_survivors_1_5_GSE54514.RData",
    "data/GO/statistics_score_metrics_non-survivors_1_5_GSE54514.RData"
  ),
  group_labels = c("HC 1-5", "S 1-5", "NS 1-5"),
  levels_order = c("HC 1-5", "S 1-5", "NS 1-5"),
  output_file = "GSE54514_GO_stat.svg",
  colors = c("yellowgreen", "snow4", "violetred")
)

# For GSE95233
plot_statistics(
  files = c(
    "data/KEGG/statistics_score_metrics_survivors_1_2_GSE95233.RData",
    "data/KEGG/statistics_score_metrics_survivors_1_3_GSE95233.RData",
    "data/KEGG/statistics_score_metrics_non-survivors_1_2_GSE95233.RData",
    "data/KEGG/statistics_score_metrics_non-survivors_1_3_GSE95233.RData"
  ),
  group_labels = c("S 1-2","S 1-3","NS 1-2","NS 1-3"),
  levels_order = c("S 1-2","S 1-3","NS 1-2","NS 1-3"),
  output_file = "GSE95233_KEGG_stat.svg",
  colors = c("snow4","snow4", "violetred", "violetred")
)

plot_statistics(
  files = c(
    "data/GO/statistics_score_metrics_survivors_1_2_GSE95233.RData",
    "data/GO/statistics_score_metrics_survivors_1_3_GSE95233.RData",
    "data/GO/statistics_score_metrics_non-survivors_1_2_GSE95233.RData",
    "data/GO/statistics_score_metrics_non-survivors_1_3_GSE95233.RData"
  ),
  group_labels = c("S 1-2","S 1-3","NS 1-2","NS 1-3"),
  levels_order = c("S 1-2","S 1-3","NS 1-2","NS 1-3"),
  output_file = "GSE95233_GO_stat.svg",
  colors = c("snow4","snow4", "violetred", "violetred")
)



# Load data KEGG
#==========================================================================================================

main_dir <- paste0(path_to_data, "/GSE54514/KEGG")
load_rdata_from_subdirectories(main_dir)

score_metrics_survivors_1_2_gse54514 <- score_metrics_survivors_1_2
score_metrics_survivors_1_3_gse54514 <- score_metrics_survivors_1_3
`score_metrics_non-survivors_1_2_gse54514` <- `score_metrics_non-survivors_1_2`
`score_metrics_non-survivors_1_3_gse54514` <- `score_metrics_non-survivors_1_3`
main_dir <- paste0(path_to_data, "/GSE48080/KEGG")
load_rdata_from_subdirectories(main_dir)

main_dir <- paste0(path_to_data, "/GSE95233/KEGG")
load_rdata_from_subdirectories(main_dir)

# Load data GO
#==========================================================================================================
main_dir <- paste0(path_to_data, "/GSE54514/GO")
load_rdata_from_subdirectories(main_dir)

score_metrics_survivors_1_2_gse54514 <- score_metrics_survivors_1_2
score_metrics_survivors_1_3_gse54514 <- score_metrics_survivors_1_3
`score_metrics_non-survivors_1_2_gse54514` <- `score_metrics_non-survivors_1_2`
`score_metrics_non-survivors_1_3_gse54514` <- `score_metrics_non-survivors_1_3`

main_dir <- paste0(path_to_data, "/GSE48080/GO")
load_rdata_from_subdirectories(main_dir)

main_dir <- paste0(path_to_data, "/GSE95233/GO")
load_rdata_from_subdirectories(main_dir)


# Figure 2B
#==========================================================================================================

# Get the data
metric_list <- list(
  "HC 1-5 (GSE54514)"           = score_metrics_control_1_5$Median,
  "S 1-2 (GSE54514)"            = score_metrics_survivors_1_2_gse54514$Median,
  "S 1-3 (GSE54514)"            = score_metrics_survivors_1_3_gse54514$Median,
  "S 1-4 (GSE54514)"            = score_metrics_survivors_1_4$Median,
  "S 1-5 (GSE54514)"            = score_metrics_survivors_1_5$Median,
  "S 1-2 (GSE95233)"            = score_metrics_survivors_1_2$Median,
  "S 1-3 (GSE95233)"            = score_metrics_survivors_1_3$Median,
  "S 1-7 (GSE48080)"            = score_metrics_survivors_1_7$Median,
  "NS 1-2 (GSE54514)"           = `score_metrics_non-survivors_1_2_gse54514`$Median,
  "NS 1-3 (GSE54514)"           = `score_metrics_non-survivors_1_3_gse54514`$Median,
  "NS 1-4 (GSE54514)"           = `score_metrics_non-survivors_1_4`$Median,
  "NS 1-5 (GSE54514)"           = `score_metrics_non-survivors_1_5`$Median,
  "NS 1-2 (GSE95233)"           = `score_metrics_non-survivors_1_2`$Median,
  "NS 1-3 (GSE95233)"           = `score_metrics_non-survivors_1_3`$Median,
  "NS 1-7 (GSE48080)"           = `score_metrics_non-survivors_1_7`$Median
)

data <- data.frame(
  Median = unlist(metric_list, use.names = FALSE),
  Group = factor(rep(names(metric_list), times = sapply(metric_list, length)),
                 levels = names(metric_list))
)


data$Category <- factor(
  ifelse(data$Group %in% c("HC 1-5 (GSE54514)"), "HC",
         ifelse(data$Group %in% c("S 1-2 (GSE54514)", "S 1-3 (GSE54514)",
                                  "S 1-4 (GSE54514)", "S 1-5 (GSE54514)",
                                  "S 1-2 (GSE95233)", "S 1-3 (GSE95233)",
                                  "S 1-7 (GSE48080)"), "S",
                "NS")),
  levels = c("HC", "S", "NS")
)

data$Group <- factor(data$Group, levels = unique(data$Group))

# Plot
raincloud_plot <- ggplot(data, aes(x = Group, y = Median, fill = Category, color = Category)) +

  # Half-violin density plot
  stat_halfeye(
    adjust = 0.5,
    justification = -0.2,
    .width = 0,
    point_colour = NA
  ) +

  # Boxplot
  geom_boxplot(
    width = 0.15,
    alpha = 0.5,
    outlier.color = NA
  ) +

  # Jittered points
  geom_jitter(
    width = 0.01,
    alpha = 0.2,
    size = 0.1
  ) +

  # Axis and Aesthetics
  scale_fill_manual(values = c("HC" = 'yellowgreen',
                               "S" =  "snow4",
                               "NS" = "violetred")) +
  scale_color_manual(values = c("HC" = 'yellowgreen',
                                "S" =  "snow4",
                                "NS" = "violetred")) +
  labs(
    x = "",
    y = "Median rewiring scores",
    fill = "Category"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11), # Rotate x-axis labels
    legend.title = element_blank(),
    axis.text.y=element_text(size = 11),
    legend.position = "bottom",
    plot.margin = margin(10, 10, 10, 50)
  )
ggsave("raincloud_plot_GO.svg", plot = raincloud_plot, device = "svg", width = 9, height = 6)


# Figure 3
#==========================================================================================================
df_list <- list(score_metrics_control_1_5, score_metrics_survivors_1_5, score_metrics_survivors_1_7,
                score_metrics_survivors_1_3, `score_metrics_non-survivors_1_5`, `score_metrics_non-survivors_1_7`,`score_metrics_non-survivors_1_3`)

groups_factor <- c('HC', 'S', 'S', 'S', 'NS', 'NS','NS')

top_n <- 10
top_pathways_list <- list()

#Get top 10 pathways per group
for (df in df_list) {
  df_sorted <- df %>% arrange(desc(Median))
  top_pathways <- head(df_sorted$Pathway, top_n)
  top_pathways_list[[length(top_pathways_list) + 1]] <- top_pathways
}

# Get union of pathways
all_pathways <- unique(unlist(top_pathways_list))

# Create heatmap
heatmap_data <- data.frame(Pathway = all_pathways)

# Get median values for all pathways
for (i in seq_along(df_list)) {
  df <- df_list[[i]]

  pathway_values <- sapply(all_pathways, function(pathway) {
    median_value <- df %>% filter(Pathway == pathway) %>% pull(Median)
    if (length(median_value) == 0) {
      return(NA)
    } else {
      return(median_value)
    }
  })

  heatmap_data <- cbind(heatmap_data, pathway_values)
}

colnames(heatmap_data)[-1] <- c('HC 1-5 (GSE54514)', 'S 1-5 (GSE54514)', 'S 1-7 (GSE48080)', 'S 1-3 (GSE95233)', 'NS 1-5 (GSE54514)', 'NS 1-7 (GSE48080)', 'NS 1-3 (GSE95233)')
heatmap_data_melted <- melt(heatmap_data, id.vars = "Pathway")

# Get colors
jet.colors <- function(n) {
  colorRampPalette(c("#86d780", "#ffea70", "coral", 'red'))(n)
}

# Plot heatmap
p <- ggplot(heatmap_data_melted, aes(x = variable, y = Pathway, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c(jet.colors(100)),
    na.value = "white",
    limits = c(min(heatmap_data_melted$value, na.rm = TRUE),
               max(heatmap_data_melted$value, na.rm = TRUE))
  ) +
  theme_minimal() +
  labs(x = "", y = "Pathways", fill = "Median score") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
        axis.text.y = element_text(size = 8))
p
ggsave("heatmap_plot_top10_GO.svg", plot = p, device = "svg", width = 8, height = 11)
library(writexl)

write_xlsx(heatmap_data_melted, "heatmap_data_GO.xlsx")

# Figure 4
#==========================================================================================================
df_list <- list(score_metrics_control_1_5, score_metrics_survivors_1_2_gse54514, score_metrics_survivors_1_3_gse54514,
                score_metrics_survivors_1_4, score_metrics_survivors_1_5, score_metrics_survivors_1_7, score_metrics_survivors_1_2,
                score_metrics_survivors_1_3,  `score_metrics_non-survivors_1_2_gse54514`,  `score_metrics_non-survivors_1_3_gse54514`,
                `score_metrics_non-survivors_1_4`, `score_metrics_non-survivors_1_5`, `score_metrics_non-survivors_1_7`, `score_metrics_non-survivors_1_2`,`score_metrics_non-survivors_1_3`)

extract_median_of_medians <- function(df, dataset, group, day_comparison) {
  median_value <- median(df$Median, na.rm = TRUE)
  data.frame(
    median_of_medians = median_value,
    dataset = dataset,
    group = group,
    day_comparison = day_comparison
  )
}
datasets <- c("GSE54514", "GSE54514", "GSE54514", "GSE54514", "GSE54514", "GSE48080", "GSE95233", "GSE95233",
              "GSE54514", "GSE54514", "GSE54514", "GSE54514", "GSE48080", "GSE95233", "GSE95233")
groups <- c("control", "survivors", "survivors", "survivors", "survivors", "survivors", "survivors", "survivors",
            "non-survivors", "non-survivors", "non-survivors", "non-survivors", "non-survivors", "non-survivors", "non-survivors")
day_comparisons <- c("1-5", "1-2", "1-3", "1-4", "1-5", "1-7", "1-2", "1-3", "1-2", "1-3", "1-4", "1-5", "1-7", "1-2", "1-3")

summary_df_go <- do.call(rbind, lapply(1:length(df_list), function(i) {
  extract_median_of_medians(df_list[[i]], datasets[i], groups[i], day_comparisons[i])
}))

shape_map <- c("GSE54514" = 22, "GSE48080" = 21, "GSE95233" = 24)
color_map <- c("control" = "yellowgreen", "survivors" = "snow4", "non-survivors" = "violetred")

# Plot
ggplot(summary_df_go, aes(x = day_comparison, y = median_of_medians)) +
  geom_point(aes(shape = dataset, fill = group), size = 4, color = "black") +
  scale_shape_manual(values = shape_map) +
  scale_fill_manual(values = color_map) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, size = 5)),
    shape = guide_legend(override.aes = list(fill = NA))
  ) +
  labs(
    x = "",
    y = "",
    fill = "Group",
    shape = "Dataset"
  ) +
  theme_minimal()
ggsave("median_of_medians_plot_GO.svg", device = "svg", width = 6, height = 4)
