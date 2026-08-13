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
datasets        <- 'Dataset_name'            # e.g. GSEXXXXX
database        <- 'KEGG'                    # KEGG or GO
output_dir      <- 'results'                 # directory to export results
groups          <- c('Healthy', 'Patients')  # groups labels to analyze
timepoint_pairs <- c(1,2)                    # A pair of points e.g. day 1 and day 2
load('ids.RData')                            # a dataframe with metadata including time_col, subject_ids, 
                                             # sample_ids, and group_col

# Get database
#==========================================================================================================
# To get the latest version of each database re-run the individual scripts for KEGG and GO export data in /R.

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
clean_matrix <- prepare_gene_intensities(
  norm = FALSE,                                   # whether to quantile normalize or not
  log = FALSE,                                    # whether to log2-transform
  series_matrix_file = 'series_matrix.txt',       # the series matrix file
  platform_file = 'platform.txt',                 # the platform annotation file
  probe_col = "ID",                               # column in platform file with probe IDs
  symbol_col = "Symbol",                          # column with gene symbols
  entrez_col = "Entrez_Gene_ID",                  # columns with Entrez IDs
  series_skip = 0,                                # header lines to skip in series matrix
  platform_skip = 0,                              # header lines to skip in platform file
  series_sep = "\t",
  platform_sep = "\t",
  trim_cols = TRUE
)

gi=as.matrix(clean_matrix[,-1])

# Run pipeline
#==========================================================================================================
for (group in groups) {
  for (tp_pair in timepoint_pairs) {
    tp <- as.numeric(tp_pair)
    labels <- as.character(tp)
    
    net <- new("gene_network",
               gene_intensities = gi,
               group_name       = group,
               timepoints       = tp,
               metadata         = ids,
               output_dir       = output_dir,
               dataset          = dataset,
               database         = database,
               pathways         = pathways,
               group_col        = "Group",
               subjects_ids     = "id",
               sample_ids       = "ID",
               time_col         = "Day",
               time_labels      = labels)
    
    run_analysis(net)
  }
}