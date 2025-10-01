# Libraries
#==========================================================================================================
suppressPackageStartupMessages({
  library(WGCNA)
  library(data.table)
  library(igraph)
  library(hash)
  library(stringi)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(moments)
  library(stats)
  library(stringr)
  library(NetworkRewiring)
})

# Select data
#==========================================================================================================
datasets <- c("GSE54514", "GSE95233", "GSE48080")
databases <- c("KEGG", "GO")

dataset = datasets[1]
database = databases[1]

output_dir = '/path/to/results'
load(paste('examples',"data", dataset, "ids.RData", sep = "/"))

# Get database
#==========================================================================================================
if(database=="KEGG"){
  load(paste('data', "KEGGpathways2024_Entrez.RData", sep = "/"))
  pathways <- pathways2023_Entrez
  rm(pathways2023_Entrez)
}else{
  load(paste('data',"biological_processes.RData", sep = "/"))
  pathways <- processes
  pathways <- subset(pathways, Total > 2 & Total <=500)
  rm(processes)
}

# Get gene expression
#==========================================================================================================
if (dataset == "GSE54514") {

  possible_groups <- c("control", "survivors", "non-survivors")
  combinations <- t(combn(c("1", "2", "3", "4", "5"), 2))

  clean_matrix <- prepare_gene_intensities(
    norm = FALSE,
    log = TRUE,
    series_matrix_file = paste('examples',"data", dataset, "GSE54514_series_matrix.txt", sep = "/"),
    platform_file = paste('examples',"data", dataset, "GPL6947-13512.txt", sep = "/"),
    probe_col = "ID",
    symbol_col = "Symbol",
    entrez_col = "Entrez_Gene_ID",
    series_skip = 71,
    platform_skip = 30,
    series_sep = "\t",
    platform_sep = "\t",
    trim_cols = TRUE
  )

} else if (dataset == "GSE95233") {

  possible_groups <- c("survivors", "non-survivors")
  combinations <- data.frame(c(1,1),c(2,3))

  clean_matrix <- prepare_gene_intensities(
    norm = FALSE,
    log = TRUE,
    series_matrix_file = paste('examples',"data", dataset,"GSE95233_series_matrix.txt", sep = "/"),
    platform_file = paste('examples',"data", dataset, "GPL570-55999.txt", sep = "/"),
    probe_col = "ID",
    symbol_col = "Gene.Symbol",
    entrez_col = "ENTREZ_GENE_ID",
    series_skip = 71,
    platform_skip = 16,
    series_sep = "\t",
    platform_sep = "\t",
    gene_synonyms = TRUE

  )
} else if (dataset == "GSE48080") {

  possible_groups <- c("survivors", "non-survivors")
  combinations <- data.frame(1,7)

  clean_matrix <- prepare_gene_intensities(
    norm = FALSE,
    log = TRUE,
    series_matrix_file = paste('examples',"data", dataset, "GSE48080_series_matrix.txt", sep = "/"),
    platform_file = paste('examples',"data", dataset, "GPL4133-12599.txt", sep = "/"),
    probe_col = "ID",
    symbol_col = "GENE_SYMBOL",
    entrez_col = "GENE",
    series_skip = 72,
    platform_skip = 22,
    series_sep = "\t",
    platform_sep = "\t",
    numeric_cols = c(1,3)
  )
}

# Prepare gene expression
#==========================================================================================================
gi=as.matrix(clean_matrix[,-1])

# Run pipeline
#==========================================================================================================
for (group in possible_groups) {

  tp_list <- if (group == "control") list(c(1, 5)) else split(combinations, seq(nrow(combinations)))

  for (tp_pair in tp_list) {
    tp <- as.numeric(tp_pair)
    labels <- as.character(tp)

    net <- new("gene_network",
               gene_intensities = gi,
               group_name = group,
               timepoints = tp,
               metadata = ids,
               output_dir = output_dir,
               dataset = dataset,
               database = database,
               pathways = pathways,
               group_col = "Group",
               subjects_ids = "id",
               sample_ids = "ID",
               time_col = "Day",
               time_labels = labels)

    run_analysis(net)
  }
}




