#' Get GO Biological Processes
#'
#' Reads a gene-to-GO mapping file and constructs a table of biological
#' processes with their associated gene sets.
#'
#' @details
#' The input file `gene_attribute_edges.txt` must contain at least the
#' following columns:
#' - `GO.Biological.Process`: Biological process name
#' - `Gene.ID`: Entrez gene ID or equivalent identifier
#'
#' The function aggregates all genes per process, counts them, and
#' outputs a `data.frame` with:
#' - `Total`: Number of genes in the process
#' - `Process`: Process name
#' - `Entrez_ID`: Comma-separated list of genes
#' - Extra columns: Each gene in a separate column (wide format)
#'
#' Finally, the result is saved as `biological_processes.RData`.
#'
#' @importFrom stringi stri_list2matrix
#' @export
get_biological_processes_GO <- function(input_file = "gene_attribute_edges.txt",
                               output_file = "biological_processes.RData") {
  # Load GO annotation edges
  go <- read.table(input_file, header = TRUE, sep = "\t", quote = "")

  # Unique biological processes
  process_names <- data.frame(Process = unique(go$GO.Biological.Process))

  # Map each process to its associated gene IDs
  biological_processes <- lapply(process_names$Process, function(proc) {
    go[go$GO.Biological.Process == proc, "Gene.ID"]
  })
  names(biological_processes) <- process_names$Process

  # Build process summary table
  totals <- sapply(biological_processes, length)
  gene_strings <- sapply(biological_processes, function(ids) {
    paste(ids, collapse = ",")
  })

  processes <- data.frame(
    Total = totals,
    Process = names(totals),
    Entrez_ID = gene_strings,
    stringsAsFactors = FALSE
  )

  # Expand gene lists into wide matrix (one column per gene)
  entrez_matrix <- stringi::stri_list2matrix(strsplit(processes$Entrez_ID, ","))
  entrez_df <- as.data.frame(t(entrez_matrix), stringsAsFactors = FALSE)

  processes <- cbind(processes, entrez_df)

  # Ensure numeric Total and sort by size
  processes$Total <- as.numeric(processes$Total)
  processes <- processes[order(processes$Total), ]

  rownames(processes) <- NULL

  # Save to file
  save(processes, file = output_file)

  message("Biological processes saved to ", output_file)
  invisible(processes)
}
