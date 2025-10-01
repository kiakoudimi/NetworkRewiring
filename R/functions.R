#' Gene co-expression network class
#'
#' A class to represent a pathway-specific gene co-expression network for estimating network rewiring in two time points.
#'
#' @slot group_name Name of the group to analyze (e.g. HC, Patients)
#' @slot timepoints Vector of timepoints to analyze (numeric or character)
#' @slot gene_intensities Matrix of gene expression data, where rows are genes and columns are sample or subject IDs
#' @slot metadata Data frame of sample metadata describing the samples, including group, time, subject IDs
#' @slot pathways Data frame of pathway/gene sets
#' @slot output_dir Path to the output directory
#' @slot dataset The name or accession ID of the dataset
#' @slot database Name of database (KEGG or GO)
#' @slot group_col Column name in `metadata` representing group
#' @slot subjects_ids Column name in `metadata` representing subjects IDs
#' @slot sample_ids Column name in `metadata` representing sample IDs
#' @slot time_col Column name in `metadata` representing timepoints
#' @slot time_labels Optional labels for timepoints to use in file names
#' @export
setClass("gene_network",
         slots = list(
           group_name = "character",
           timepoints = "vector",
           gene_intensities = "matrix",
           metadata = "data.frame",
           pathways = "data.frame",
           output_dir = "character",
           dataset = "character",
           database = "character",
           group_col = "character",
           subjects_ids = "character",
           sample_ids = "character",
           time_col = "character",
           time_labels = "vector"
         ))


#' Generic for select_group_intensities
#'
#' @param object An object
#' @export
setGeneric("select_group_intensities", function(object,...) {
  standardGeneric("select_group_intensities")
})

#' Generic for network_construction
#'
#' @param object An object
#' @export
setGeneric("network_construction", function(object) {
  standardGeneric("network_construction")
})

#' Generic for compute_rewiring
#'
#' @param object An object
#' @export
setGeneric("compute_rewiring", function(object,...) {
  standardGeneric("compute_rewiring")
})

#' Select group intensities
#'
#' Extract gene intensities for a given group across two timepoints.
#'
#' @param object A \code{gene_network} object.
#'
#' @return A list with two matrices: gene intensities for time_1 and time_2.
#' @export
setMethod("select_group_intensities", "gene_network",
          function(object) {
            group_name <- object@group_name
            timepoints <- object@timepoints
            ids <- object@metadata
            gene_intensities <- object@gene_intensities

            # Required metadata columns
            required_cols <- c(object@group_col, object@time_col, object@subjects_ids, object@sample_ids)
            missing_cols <- setdiff(required_cols, colnames(ids))
            if (length(missing_cols) > 0) {
              stop(sprintf("Metadata is missing required column(s): %s",
                           paste(missing_cols, collapse = ", ")))
            }
            #Ensure at least two timepoints are given
            if (length(timepoints) < 2) stop("At least two timepoints are required")

            # Labels for output
            labels <- if (!is.null(object@time_labels)) object@time_labels else timepoints

            # Extract subjects per timepoint
            subjects_by_time <- lapply(timepoints, function(tp) {
              subset <- ids[
                ids[[object@group_col]] == group_name &
                  ids[[object@time_col]] == tp, , drop = FALSE
              ]
              if (nrow(subset) == 0) {
                warning(sprintf("No samples found for group '%s' at timepoint '%s'", group_name, tp))
                return(list(subjects = character(0), samples = character(0)))
              }
              list(
                subjects = subset[[object@subjects_ids]],         # subject IDs
                samples  = subset[[object@sample_ids]]    # acquisition/sample IDs
              )
            })
            names(subjects_by_time) <- labels

            # keep only subjects present at all timepoints
            common_subjects <- Reduce(intersect, lapply(subjects_by_time, `[[`, "subjects"))
            subjects_by_time <- lapply(subjects_by_time, function(x) {
              keep_idx <- which(x$subjects %in% common_subjects)
              list(subjects = x$subjects[keep_idx],
                   samples  = x$samples[keep_idx])
            })

            intensities_by_time <- lapply(subjects_by_time, function(x) {
              if (length(x$samples) == 0) return(NULL)
              sample_ids <- intersect(x$samples, colnames(gene_intensities))
              if (length(sample_ids) == 0) {
                warning("No matching sample IDs found in gene_intensities for this timepoint")
                return(NULL)
              }

              mat <- gene_intensities[, sample_ids, drop = FALSE]

              # Convert to numeric while keeping rownames
              mat <- matrix(as.numeric(mat), nrow = nrow(mat), ncol = ncol(mat),
                            dimnames = list(rownames(mat), colnames(mat)))

              return(mat)
            })
          })

#' Network construction
#'
#' Constructs adjacency networks for a given group and two timepoints.
#'
#' @param object A \code{gene_network} object.
#' @export
setMethod("network_construction", "gene_network",
          function(object) {
            group <- object@group_name
            timepoints <- object@timepoints
            output_dir <- object@output_dir
            dataset <- object@dataset
            pathways <- object@pathways
            database <- object@database

            #Ensure at least two timepoints are given
            if (length(timepoints) < 2) stop("At least two timepoints are required")

            # Labels for output
            labels <- if (!is.null(object@time_labels)) object@time_labels else timepoints
            tp1_label <- labels[1]
            tp2_label <- labels[2]

            # Get expression matrices for the group at each timepoint
            data_list <- select_group_intensities(object)
            group_1 <- data_list[[tp1_label]]
            group_2 <- data_list[[tp2_label]]

            n_pathways <- nrow(pathways)

            for (i in seq_len(n_pathways)) {
              cat(sprintf("\rProcessing pathway %d/%d ...", i, n_pathways))

              # Extract pathway genes
              if (database == "KEGG") {
                path_name <- pathways$Pathway[i]
                path <- pathways[pathways$Pathway == path_name, -c(1:2)]

              } else {
                path_name <- pathways$Process[i]
                path <- pathways[pathways$Process == path_name, -c(1:3)]
              }

              path <- path[, colSums(is.na(path)) < 1]
              path_genes <- unlist(path)
              path_genes <- intersect(path_genes, rownames(group_1))

              # Subset expression matrices
              group_1_path <- group_1[path_genes, , drop = FALSE]
              group_2_path <- group_2[path_genes, , drop = FALSE]

              # Ensure at least two genes are present
              if (nrow(group_1_path) >= 2 && nrow(group_2_path) >= 2) {

                # Compute adjacency
                apath1 <- adjacency(t(group_1_path), type = "signed", power = 1)
                apath2 <- adjacency(t(group_2_path), type = "signed", power = 1)

                # Clean path name for filenames
                if (grepl("/", path_name)){
                  path_name <- gsub("/", ",", path_name)
                }

                # Create output directory
                base_dir <- file.path(output_dir, dataset, database,
                                      paste0(group, "_", tp1_label, "_", tp2_label), "edge_lists")
                if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)

                edgefile1 <- file.path(base_dir, paste0(i, ".1.", group, "_", tp1_label, "_", path_name, ".txt"))
                edgefile2 <- file.path(base_dir, paste0(i, ".2.", group, "_", tp2_label, "_", path_name, ".txt"))

                # Export networks
                exportNetworkToCytoscape(apath1, edgeFile = edgefile1, weighted = TRUE, threshold = 0,
                                         nodeNames = colnames(apath1))
                exportNetworkToCytoscape(apath2, edgeFile = edgefile2, weighted = TRUE, threshold = 0,
                                         nodeNames = colnames(apath2))

              }
            }

            cat("\nNetwork construction completed.\n")
          })

#' Network rewiring
#'
#' Computes network rewiring for a given group and two timepoints.
#'
#' @param object A \code{gene_network} object.
#' @param edgelist_path Path to the directory containing edge list files
#' @param filenames Filenames refering to the pathway names
#' @return Export the rewiring scores for each pathway.
#' @export
setMethod("compute_rewiring", "gene_network",
          function(object, edgelist_path, filenames){
            group <- object@group_name
            output_dir <- object@output_dir
            database <- object@database
            timepoints <- object@timepoints

            # Labels for output
            labels <- if (!is.null(object@time_labels)) object@time_labels else timepoints
            tp1_label <- labels[1]
            tp2_label <- labels[2]

            output.folder <- file.path(output_dir, dataset, database,
                                  paste0(group, "_", tp1_label, "_", tp2_label), "dynet_score")

            if (!dir.exists(output.folder)) {
              dir.create(output.folder, recursive = TRUE)
            }

            k <- 1
            n_iters <- length(seq(1, length(filenames), by = 2))
            for(p in seq(1,length(filenames), by = 2)){
              cat(paste0("\rComputing network rewiring ", k, "/", length(filenames)/2, " ..."), end='\r')
              total_edge_lists <- lapply(paste(edgelist_path, filenames, sep="/")[p:(p+1)], read.table, sep="\t",stringsAsFactors = FALSE, header=T)
              names(total_edge_lists)<-c(letters[seq(from = 1, to = length(total_edge_lists))])

              if(nrow(total_edge_lists$a)!=0){

                #Set cutoff (for edges)
                cutoff=0

                top_mes<-NULL

                for (i in 1:length(total_edge_lists)){
                  total_edge_lists[[i]]<-total_edge_lists[[i]][,1:3]
                  total_edge_lists[[i]]$edge<-paste(total_edge_lists[[i]]$fromNode, total_edge_lists[[i]]$toNode, sep="_")
                  top_mes<-union(top_mes, total_edge_lists[[i]][total_edge_lists[[i]]$weight > cutoff,]$edge)
                }

                #Calculation of rewiring scores
                dn_scores<-do.call(rbind.data.frame,compare_multiple_networks(total_edge_lists,top_mes))

                #---------------------------------------------------------------------------
                #Export Data
                filename<-paste0(output.folder,"/", sub(paste0("1.",group,"_day_1*."),"",filenames[p]) )
                write.table(dn_scores,file=filename, sep="\t", quote = F, row.names=FALSE)

              }
              k=k+1
            }
            cat("\nNetwork rewiring completed.\n")
          })

#' Network rewiring function
#' @param x The edge lists as obtained from \code{compute_rewiring}
#' @param top_mes Union lists as obtained from \code{compute_rewiring}
#' @return A numeric matrix of the same dimensions with or without normalization/log-transform.
#' @export
compare_multiple_networks <- function(x,top_mes) {

  # Find all nodes in networks
  megenes<-unique(sort(unlist(strsplit(top_mes,"_"))))
  allAdg<-NULL
  for (net in names(x)){

    x_alt<-x[[net]][x[[net]]$edge %in% top_mes,]

    #edge list to adjacency matrix
    b<-data.table(from=x_alt$fromNode,to=x_alt$toNode,score=x_alt$weight)

    apnodes<-setdiff(megenes,unique(sort(c(b$from,b$to))))

    g=graph_from_data_frame(b,directed=FALSE)
    e=as_adjacency_matrix(g,type="both",names=TRUE,sparse=FALSE,attr="score");
    if(length(apnodes)!=0){
      genes<-rownames(e)
      for(i in 1:length(apnodes)){
        e<-rbind(e,rep(0,dim(e)[2]))
      }
      for(i in 1:length(apnodes)){
        e<-cbind(e,rep(0,dim(e)[1]))
      }

      colnames(e)<-c(genes,apnodes)
      rownames(e)<-c(genes,apnodes)
    }
    e<-e[order(rownames(e)), ]
    e<-e[, order(colnames(e))]
    e[e == ""]<-0

    for(i in 1:dim(e)[1]){
      allAdg[[net]][[rownames(e)[i]]]<-as.numeric(e[i,])
    }

  }

  #find mean of non zero weights over all states
  meantable<-NULL
  for (g in megenes){
    tmpar<-NULL
    for (net in names(allAdg)){
      tmpar<-rbind(tmpar,allAdg[[net]][[g]])
    }
    is.na(tmpar) <- tmpar==0
    meantable[[g]]<-colMeans(tmpar, na.rm=TRUE)
    meantable[[g]][is.na(meantable[[g]])] <- 0
  }

  #standartized weights divided by mean
  for (net in names(allAdg)){
    for (g in megenes){
      allAdg[[net]][[g]]<-allAdg[[net]][[g]]/meantable[[g]]
      allAdg[[net]][[g]][is.na(allAdg[[net]][[g]])] <- 0
    }
  }

  #compute centroid matrix
  for (g in megenes){
    tmpar<-NULL
    for (net in names(allAdg)){
      tmpar<-rbind(tmpar,allAdg[[net]][[g]])
    }
    # print(tmpar)
    meantable[[g]]<-colMeans(tmpar, na.rm=TRUE)
  }

  #eukleidian distance to centriod
  dn_ar<-NULL
  for (net in names(allAdg)){
    dn<-NULL
    for (g in megenes){
      dn<-c(dn,sum((allAdg[[net]][[g]]-meantable[[g]])^2))
    }
    dn_ar<-cbind(dn_ar,dn)
  }
  colnames(dn_ar)<-names(allAdg)

  dn_ar1<-dn_ar
  dn_scores<-list()
  edges<-data.table(do.call(rbind,strsplit(top_mes,"_")))
  colnames(edges)<-c("from","to")
  for (i in 1:length(megenes)){
    dn_scores[[megenes[i]]]$gene <- megenes[i]
    dn_scores[[megenes[i]]]$edges<-
      dim(edges[(edges$from == megenes[i]) | (edges$to == megenes[i]),])[1]
    # MZ: added cytoscape version
    dn_scores[[megenes[i]]]$score_cyto<-sum(dn_ar1[i,])/(dim(dn_ar)[2]-1)
    dn_scores[[megenes[i]]]$score_cor_cyto<-
      dn_scores[[megenes[i]]]$score_cyto/dn_scores[[megenes[i]]]$edges
  }
  return(dn_scores)
}

# Check if the data needs log-transformation
#' @param intensities A numeric matrix of gene intensities.
#' @return TRUE or FALSE of whether to proceed with log-transform
#' @export
checkLog <- function(intensities){
  toLog <- FALSE
  if (max(intensities, na.rm = TRUE) > 25) toLog <- TRUE
  return(toLog)
}

# Normalize and/or log-transform expression data
#' @param intensities A numeric matrix of gene intensities.
#' @param norm Whether to normalize the data or not.
#' @param log Whether to log-transform the data or not.
#' @return A numeric matrix of the same dimensions with ot without normalization/log-transform.
#' @export
handleData <- function(intensities, norm = FALSE, log = FALSE){
  if (norm && log){
    cat("Normalizing and doing a log2 on the data.\n")
    intensities <- normalizeQuantiles(intensities)
    intensities <- log2(intensities)
  } else if (!norm && log){
    cat("Doing a log2 only on the data.\n")
    intensities <- log2(intensities)
  } else if (norm && !log){
    cat("Un-log, normalize, and re-log the data.\n")
    intensities <- 2^intensities
    intensities <- normalizeQuantiles(intensities)
    intensities <- log2(intensities)
  }
  return(intensities)
}

# Remove duplicate genes (general)
#' @param eset A numeric matrix of gene intensities.
#' @param column_of_symbol The columns name of gene symbols
#' @return A numeric matrix of the same dimensions without duplicate acquisitions.
#' @export
remove_duplicate_genes <- function(eset, column_of_symbol = "GeneSymbol", method = "mean") {
  symbols <- eset[[column_of_symbol]]
  unique_symbols <- unique(symbols)

  result <- lapply(unique_symbols, function(sym) {
    rows <- which(symbols == sym)
    if (length(rows) == 1) {
      return(eset[rows, , drop = FALSE])
    } else {
      exprs <- eset[rows, -which(names(eset) == column_of_symbol), drop = FALSE]
      if (method == "mean") {
        exprs_mean <- colMeans(exprs)
      }
      return(cbind(GeneSymbol = sym, t(exprs_mean)))
    }
  })

  result <- do.call(rbind, result)
  rownames(result) <- result[, column_of_symbol]
  result <- as.data.frame(result, stringsAsFactors = FALSE)

  return(result)
}

# Handle negatives
#' @param intensities A numeric matrix of gene intensities.
#' @return A numeric matrix or vector of the same dimensions with all values shifted to be positive.
#' @export
cureNegative <- function(intensities){
  cat("Handling negatives on the data.\n")
  min <- abs(min(intensities)) + 0.001
  intensities <- intensities + min
  return(intensities)
}

#' Prepare gene intensities
#'
#' Cleans a raw gene expression matrix, handles multiple probes mapping to the same gene, ensures data are normalized and log-transformed, and maps gene symbols to Entrez IDs.
#'
#' @param norm Logical. If TRUE, the expression data will be quantile-normalized. Default is FALSE.
#' @param log Logical. If TRUE, the method will log2-transform the data if necessary. Default FALSE
#' @param series_matrix_file Character. Path to the raw series matrix file. Optional if the matrix is already in `gene_intensities`.
#' @param platform_file Character. Path to the platform file mapping probes to gene symbols and Entrez IDs. Required for mapping.
#' @param probe_col Character. Column name in the platform file containing probe IDs. Default is "ID".
#' @param symbol_col Character. Column name in the platform file containing gene symbols. Default is "Gene.Symbol".
#' @param entrez_col Character. Column name in the platform file containing Entrez IDs. Default is "ENTREZ_GENE_ID".
#' @param series_skip Integer. Number of lines to skip at the start of the series matrix file.
#' @param platform_skip Integer. Number of lines to skip at the start of the platform file.
#' @param series_sep Character. Field separator in the series matrix file. Default is tab ("\t").
#' @param platform_sep Character. Field separator in the platform file. Default is tab ("\t").
#' @param trim_cols Character vector. Column names in the platform to trim whitespace from. Default is NULL.
#' @param numeric_cols Character vector. Column names in the platform to convert to numeric. Default is NULL.
#' @param gene_synonyms Logical. If TRUE, the synonyms are handled. Default FALSE
#' @return Data frame of cleaned gene intensities with Entrez IDs as rownames
#' @export
prepare_gene_intensities <- function(norm = FALSE,
         log = FALSE,
         series_matrix_file = NULL,
         platform_file = NULL,
         probe_col = "ID",
         symbol_col = "Gene.Symbol",
         entrez_col = "ENTREZ_GENE_ID",
         series_skip = NULL,
         platform_skip = NULL,
         series_sep = "\t",
         platform_sep = "\t",
         trim_cols = FALSE,
         numeric_cols = NULL,
         gene_synonyms = FALSE) {

  library(AnnotationDbi)
  library(org.Hs.eg.db)

  # Load series matrix
  if (!is.null(series_matrix_file)) {
    series_matrix0 <- as.matrix(
      read.delim(series_matrix_file,
                 header = TRUE,
                 skip = series_skip,
                 sep = series_sep,
                 skipNul = TRUE,
                 fill = TRUE)
    )
    series_matrix0 <- series_matrix0[1:(nrow(series_matrix0)-1), ]
  }

  # Load platform
  if (!is.null(platform_file)) {
    platform <- read.delim(platform_file,
                           header = TRUE,
                           skip = platform_skip,
                           sep = platform_sep,
                           skipNul = TRUE,
                           fill = TRUE)
    platform <- data.frame(
      probeIDs = platform[[probe_col]],
      GeneSymbol = platform[[symbol_col]],
      EntrezID = platform[[entrez_col]],
      stringsAsFactors = FALSE
    )

    # Optional: convert numeric columns
    if (!is.null(numeric_cols)) {
      platform[, numeric_cols] <- lapply(platform[, numeric_cols], as.numeric)
    }

    # Optional: trim whitespace
    if (trim_cols) {
      platform <- data.frame(lapply(platform, trimws), stringsAsFactors = FALSE)
    }
  } else {
    stop("platform_file must be provided if you want to map probes to genes")
  }

  if (gene_synonyms) {
    data <- data.frame(platform, series_matrix0)
    new_records_entrez <- c()
    new_records_gene <- c()
    new_records_intens <- NULL

    for (i in seq_len(nrow(data))) {
      if (grepl("///", data$EntrezID[i])) {
        new_records_entrez <- c(new_records_entrez,
                                sub("\\///.*", "", data$EntrezID[i]),
                                sub(".*/// ", "", data$EntrezID[i]))
        new_records_gene <- c(new_records_gene,
                              sub("\\///.*", "", data$GeneSymbol[i]),
                              sub(".*/// ", "", data$GeneSymbol[i]))
        new_records_intens <- rbind(new_records_intens,
                                    c(data[i, -(2:3)]),
                                    c(data[i, -(2:3)]))
      }
    }

      if (!is.null(new_records_intens)) {
      data_upd <- data.frame(
        GeneSymbol = new_records_gene,
        EntrezID   = new_records_entrez,
        new_records_intens,
        stringsAsFactors = FALSE
      )

      remove_entry <- which(grepl("///", data$EntrezID))
      data_new <- data[-remove_entry, ]

      series_matrix <- rbind(data_new[, -1], data_upd[, -3])
    } else {
      series_matrix <- data[, -1]
    }

    platform_upd <- merge(platform, series_matrix[, c(1,2,3)], by.x = "probeIDs", by.y = "ID_REF")
    platform_upd <- platform_upd[, c(1,4,5)]
    colnames(platform_upd) <- c("probeIDs", "GeneSymbol", "EntrezID")

    platform <- platform_upd
    series_matrix0 <- series_matrix[,-c(1,2)]
  }

  temp <- as.matrix(series_matrix0[,-1])

  intensities <- matrix(as.double(temp), nrow = nrow(temp), ncol = ncol(temp))
  colnames(intensities) <- colnames(temp)
  rownames(intensities) <- series_matrix0[,1]
  remove(temp)

  # Handle normalization / log /negatives
  if (min(intensities) < 0) intensities <- cureNegative(intensities)
  if (log) toLog <- checkLog(intensities) else toLog <- FALSE
  intensities <- handleData(intensities, norm, toLog)

  eset=data.frame(probeIDs=rownames(intensities, do.NULL = FALSE, prefix=""), intensities)

  eset= merge(platform, eset, by.x = 'probeIDs', by.y = 'probeIDs')

  eset = eset[,-c(1,3)]

  gene_intensities_clean = remove_duplicate_genes(eset = eset,
                                                  column_of_symbol = 'GeneSymbol',method = 'mean')
  gene_intensities_clean <- gene_intensities_clean[-c(which(rownames(gene_intensities_clean)=="")),]

  symbols <- data.frame(GeneSymbol=rownames(gene_intensities_clean))
  entrez <- data.frame(EntrezIDs = mapIds(org.Hs.eg.db, symbols$GeneSymbol, 'ENTREZID', 'SYMBOL'))

  gene_intensities_clean <- data.frame(entrez$EntrezIDs, gene_intensities_clean)
  gene_intensities_clean <- gene_intensities_clean[complete.cases(gene_intensities_clean), ]
  rownames(gene_intensities_clean) <- gene_intensities_clean$entrez.EntrezIDs
  gene_intensities_clean <- gene_intensities_clean[,-1]
  return(gene_intensities_clean)

}

#' Compute summary statistics from rewiring scores
#'
#' @param scores The rewiring scores
#' @param filename The filename corresponding to the group and timepoint label).
#' @param dataset The dataset name
#' @param output.folder Path to output directory
#' @param pathways The pathway names
#'
#' @return Exports the statistical summary for the given dataset and group \code{output.folder}.
#'
#' @export
compute_statistics <- function(scores, filename, dataset, output.folder, pathways){
  score_metrics<-data.frame()
  path_names<-data.frame()
  sum<-data.frame()
  length<-data.frame()
  skew<-data.frame()
  top_scores<-data.frame()
  kurt<-data.frame()
  cv<-data.frame()
  for (s in 1:length(scores)){
    path<-data.frame(scores[s])
    name<-pathways[s]

    #Sum of Scores
    sum<-rbind(sum,sum(path$score_cor_cyto))

    #Total Genes
    length<-rbind(length, nrow(path))

    #Skewness
    skew<-rbind(skew, skewness(path$score_cor_cyto))

    #Kurtosis
    kurt<-rbind(kurt, kurtosis(path$score_cor_cyto))

    #Summary Metrics
    metrics<-summary(path$score_cor_cyto)
    score_metrics<-rbind(score_metrics,metrics)

    #CV
    cv<-rbind(cv,sd(path$score_cor_cyto)/mean(path$score_cor_cyto))

    #Top 5% of Scores
    top<-which(path$score_cor_cyto>quantile(path$score_cor_cyto, probs = c(0.95)))
    top<-path[c(top), c("gene","score_cor_cyto")]
    gtop<-paste(top$gene, collapse = ",")
    stop<-paste(top$score_cor_cyto, collapse = ",")
    top<-cbind(gtop,stop)
    top_scores<-rbind(top_scores,top)

    #Path Names
    path_names<-rbind(path_names,name)
  }

  score_metrics_table<-cbind(path_names,score_metrics, skew, kurt, sum, cv,length, top_scores)
  colnames(score_metrics_table)<-c("Pathway","Min.", "1st Qu.", "Median", "Mean", "3rd Qu.",
                                   "Max.", "Skewness","Kurtosis","Sum", "CV","Total Genes","Top 5% Genes", "Top 5% Score")
  assign(paste0("score_metrics_", filename), score_metrics_table)
  export_stat_path <- paste(output.folder, paste0('statistics', "_score_metrics_", filename, '_', dataset, '.RData'), sep = "/")

  # Extract the directory path
  directory_path <- dirname(export_stat_path)

  # Check if the directory exists and create it if necessary
  if (!dir.exists(directory_path)) {
    dir.create(directory_path, recursive = TRUE)
  }
  save(list = paste0("score_metrics_", filename), file = export_stat_path)
  rm(directory_path, export_stat_path)
  cat("Compute statistics completed.\n\n")
}

#' Run the full network analysis pipeline
#'
#' @param net A \code{gene_network} object containing the dataset, metadata, pathways, and configuration settings.
#'
#' @details
#' The pipeline performs the following steps:
#' \enumerate{
#'   \item Extracts settings (dataset, database, group, timepoints, labels).
#'   \item Creates output directories if they do not exist.
#'   \item Constructs networks for the specified group and timepoints via
#'     \code{\link{network_construction}}.
#'   \item Loads the exported edge list files.
#'   \item Computes rewiring scores using \code{\link{compute_rewiring}}.
#'   \item Loads rewiring scores and computes summary statistics with
#'     \code{\link{compute_statistics}}.
#' }
#'
#' The results are written into the \code{output_dir} specified in the \code{gene_network} object, under dataset- and database-specific folders.
#'
#' @return The function generates output files (edge lists, rewiring scores, and statistics) in the configured \code{output_dir}.
#' @export
run_analysis <- function(net) {

  dataset <- net@dataset
  database <- net@database
  group <- net@group_name
  time_labels <- net@time_labels
  output_dir <- net@output_dir
  name <- paste(group, time_labels[1], time_labels[2], sep = "_")
  tp <- net@timepoints

  message(sprintf(
    "Running analysis for dataset: %s | database: %s | group: %s | timepoints: %s → %s",
    dataset,
    database,
    group,
    tp[1],
    tp[2]
  ))
  flush.console()

  # Create subfolders
  edge_dir <- file.path(output_dir, dataset, database, name, "edge_lists")
  dynet_dir <- file.path(output_dir, dataset, database, name, "dynet_score")

  # Construct networks
  network_construction(net)

  # Load edge list files
  filenames <- list.files(edge_dir, pattern = "\\.txt$", full.names = FALSE)
  filenames <- str_sort(filenames, numeric = TRUE)

  # Compute rewiring scores
  compute_rewiring(net, edge_dir, filenames)

  # Load scores
  score_files <- list.files(dynet_dir, pattern = "\\.txt$", full.names = TRUE)
  score_files <- str_sort(score_files, numeric = TRUE)
  scores <- lapply(score_files, read.table, sep = "\t", stringsAsFactors = FALSE, header = TRUE)

  # Get pathway names
  if (database == "KEGG") {
    pathways <- gsub(".*?_\\d+_(.*?) - Homo sapiens \\(human\\)\\.txt", "\\1", basename(score_files))
  } else if (database == "GO") {
    pathways <- gsub(".*?_\\d+_(.*?)\\.txt", "\\1", basename(score_files))
  }

  # Compute statistics
  compute_statistics(scores, name, dataset, file.path(output_dir, dataset, database, name), pathways)
}
