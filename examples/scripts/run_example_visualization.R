# Libraries
#==========================================================================================================
library(ggplot2)
library(ggstatsplot)
library(dplyr)
library(ggdist)
library(reshape2)

library(clusterProfiler)
library(enrichplot)
library(readxl)

library(purrr)
library(stringr)
library(tidyr)
library(org.Hs.eg.db)
library(AnnotationDbi)

#Functions
#==========================================================================================================

# Get pathway statistics
load_rdata_from_subdirectories <- function(main_dir) {

  sub_dirs <- list.dirs(main_dir, full.names = TRUE, recursive = FALSE)

  for (sub_dir in sub_dirs) {
    rdata_files <- list.files(sub_dir, pattern = "\\.RData$", full.names = TRUE)

    for (file in rdata_files) {
      load(file, envir = .GlobalEnv)
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
path_to_data='...'
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
  colorRampPalette(c("#fff2a8","#d95f02","#b24a6b", "deeppink4", "#2b0033"))(n)
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

p <- ggplot(heatmap_data_melted,
            aes(x = variable, y = Pathway)) +
  geom_point(
    aes(size = value, color = value),
    shape = 16, alpha = 0.85
  ) +
  scale_color_gradientn(
    colors = jet.colors(100),
    na.value = "white"
  ) +
  scale_size(
    range = c(1, 7),   # controls bubble size
    guide = "none"     # hide size legend if redundant
  ) +
  theme_minimal() +
  labs(x = "", y = "Pathways", color = "Median score") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
    axis.text.y = element_text(size = 8)
  )

p


ggsave("heatmap_plot_top10_GO.svg", plot = p, device = "svg", width = 9.5, height = 11)
library(writexl)

write_xlsx(heatmap_data_melted, "heatmap_data_GO_21_03_2026.xlsx")

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



# Gene-level perturbations
#==========================================================================================================
ns_list <- mget(ls(pattern = "score_metrics_non-survivors_1_(5|3|7)$"))
s_list  <- mget(ls(pattern = "score_metrics_survivors_1_(5|3|7)$"))

ns_genes <- unlist(lapply(ns_list, function(df) {unlist(strsplit(df$`Top 5% Genes`, ","))}))
s_genes <- unlist(lapply(s_list, function(df) {unlist(strsplit(df$`Top 5% Genes`, ","))}))

ns_table <- sort(table(ns_genes), decreasing = TRUE)
s_table  <- sort(table(s_genes), decreasing = TRUE)

all_genes <- union(names(ns_table), names(s_table))

combined <- data.frame(
  Entrez_ID = all_genes,
  NS_freq = as.numeric(ns_table[all_genes]),
  S_freq  = as.numeric(s_table[all_genes])
)

combined$Symbol <- mapIds(
  org.Hs.eg.db,
  keys = combined$Entrez_ID,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)
combined <- combined[, c("Entrez_ID", "Symbol", "NS_freq", "S_freq")]

write_xlsx(combined, "pathway_freq_top5%_genes_KEGG_24_03_2026.xlsx")


# Pathway enrichment analysis
#==========================================================================================================
#KEGG
gene_freq <- read_excel("results/pathway_freq_top5%_genes_KEGG_24_03_2026.xlsx")
load("NetworkRewiring/data/KEGGpathways2024_Entrez.RData")
pathways <- pathways2023_Entrez
rm(pathways2023_Entrez)

term2gene <- pathways[,-2] %>%
  mutate(Pathway = sub(" - Homo sapiens \\(human\\)$", "", .[[1]])) %>%
  pivot_longer(cols = -Pathway, values_to = "gene") %>%
  filter(!is.na(gene)) %>%
  select(Pathway, gene)

#GO
gene_freq <- read_excel("results/pathway_freq_top5%_genes_GO_24_03_2026.xlsx")
load("NetworkRewiring/data/biological_processes.RData")
pathways <- processes
pathways <- subset(pathways, Total > 2 & Total <=500)
rm(processes)

term2gene <- pathways[,-c(1,3)] %>%
  pivot_longer(cols = -Process, values_to = "gene") %>%
  filter(!is.na(gene)) %>%
  select(Process, gene)

#Run OVA
percentile <- 0.9

ns_thresh <- quantile(gene_freq$NS_freq, percentile, na.rm = TRUE)
s_thresh <- quantile(gene_freq$S_freq, percentile, na.rm = TRUE)

ns_genes <- gene_freq$Entrez_ID[gene_freq$NS_freq >= ns_thresh]
s_genes <- gene_freq$Entrez_ID[gene_freq$S_freq >= s_thresh]

common_genes <- intersect(ns_genes, s_genes)
unique_ns <- setdiff(ns_genes, s_genes)
unique_s  <- setdiff(s_genes, ns_genes)

res_ns <- enricher(unique_ns, TERM2GENE = term2gene)
res_s  <- enricher(unique_s,  TERM2GENE = term2gene)
res_common <- enricher(common_genes, TERM2GENE = term2gene)

library(patchwork)

p1 <- dotplot(res_s, showCategory = 10) + ggtitle('Survivors')+
  set_enrichplot_color(type='fill', transform='log10', colors=c("#2b0033","deeppink4", "#b24a6b", 'bisque'))
p2 <- dotplot(res_ns, showCategory = 10)+ ggtitle('Non-Survivors')+
  set_enrichplot_color(type='fill', transform='log10', colors=c("#2b0033","deeppink4", "#b24a6b", 'bisque'))
p3 <- dotplot(res_common, showCategory = 10)+ ggtitle('Shared')+
  set_enrichplot_color(type='fill', transform='log10', colors=c("#2b0033","deeppink4", "#b24a6b", 'bisque'))

(p1 + p2+ p3) + plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom",axis.title.y = element_blank())
ggsave("dotplot_KEGG.svg", scale=1.3,width = 14, height = 5, device = "svg")

# Enrichment map
res_common2 <- pairwise_termsim(res_common)
res_s2 <- pairwise_termsim(res_s)
res_ns2 <- pairwise_termsim(res_ns)

ssplot(res_common2, nCluster = 5,node_label = "category", min_edge = 0.3)  +
  theme_void() +theme(legend.position = "bottom")+guides(size = "none",group = "none", fill='none' )+
  scale_x_continuous(expand = expansion(mult = 0.4)) +
  scale_y_continuous(expand = expansion(mult = 0.4))+
  set_enrichplot_color(type='color', transform='log10', colors=c("#2b0033","deeppink4", "#b24a6b", 'bisque'))

ggsave("ssplot_shared_KEGG.svg", scale=0.9,width = 15, height = 14, device = "svg")


# Heatplot
geneList_df_s <- gene_freq[gene_freq[[1]] %in% s_genes, ]
geneList_s<-geneList_df_s$S_freq
names(geneList_s) <- geneList_df_s[[1]]

geneList_df_ns <- gene_freq[gene_freq[[1]] %in% ns_genes, ]
geneList_ns<-geneList_df_ns$NS_freq
names(geneList_ns) <- geneList_df_ns[[1]]
rm(geneList_df_ns, geneList_df_s)

res_s2_symbol <- setReadable(res_s2, 'org.Hs.eg.db', 'ENTREZID')
res_ns2_symbol <- setReadable(res_ns2, 'org.Hs.eg.db', 'ENTREZID')
p1 <-heatplot(res_s2_symbol, foldChange=geneList_s, showCategory=5)
p2 <-heatplot(res_ns2_symbol, foldChange=geneList_ns, showCategory=5)
p1 + p2 + plot_layout(ncol = 1)
ggsave("heatplot_s_ns_KEGG.svg", scale=1.05,width = 23, height = 5, device = "svg")
