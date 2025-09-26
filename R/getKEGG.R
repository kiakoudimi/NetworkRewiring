# Libraries
library(magrittr)
library(dplyr)
library(org.Hs.eg.db)
library(KEGGREST)
library(stringi)

#' Get KEGG Pathways with Entrez IDs
#'
#' @param organism Character, KEGG organism code (default "hsa" for human)
#' @param output_dir Character, folder where to save the .RData file (default current working directory)
#' @param year Character or numeric, optional year for filename (default current year)
#' @return The pathways data frame is returned invisibly, also saved as .RData
#' @export
get_pathways_KEGG <- function(output_file = "KEGGpathways_Entrez.RData") {

  # Map pathways to Entrez IDs
  hsa_path_eg <- keggLink("pathway", "hsa") %>%
    tibble(pathway = ., eg = sub("hsa:", "", names(.)))

  # Pathway names/descriptions
  hsa_pathways <- keggList("pathway", "hsa") %>%
    tibble(pathway = names(.), description = .)

  # Initialize result
  pathways_Entrez <- data.frame()

  for (path in hsa_pathways$pathway) {
    entrez_ids <- hsa_path_eg[hsa_path_eg$pathway == paste0("path:", path), 2]
    entrez_ids <- paste(entrez_ids$eg, collapse = ",")
    name <- hsa_pathways[hsa_pathways$pathway == path, 2]$description
    tmp <- c(name, entrez_ids)
    pathways_Entrez <- rbind(pathways_Entrez, tmp)
  }

  colnames(pathways_Entrez) <- c("Pathway", "Entrez_ID")

  # Split Entrez_ID into separate columns
  entrez <- as.data.frame(
    t(stri_list2matrix(strsplit(pathways_Entrez$Entrez_ID, ","))),
    stringsAsFactors = FALSE
  )

  pathways_Entrez <- cbind(pathways_Entrez, entrez)

  # Clean up and save
  rm(entrez, tmp, path, name, hsa_path_eg, hsa_pathways)
  save(pathways_Entrez, file = output_file)

  cat("KEGG pathways saved to", output_file, "\n")
}
