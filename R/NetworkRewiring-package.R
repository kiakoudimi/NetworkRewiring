#'NetworkRewiring: Quantifying pathway perturbations through gene network rewiring
#'
#'A pipeline to build pathway specific gene co-expression networks from gene expression
#'data acquired in two time points and quantify how much each biological process or pathway
#'is perturbed between in those time points.
#'
#' @section Pipeline overview:
#' \enumerate{
#'  \item **Prepare gene expression data** with \code{\link{prepare_gene_intensities}}
#'    (clean proves, map gene names, normalization/log2-transform) -- if input data are already prepared skip or use the following \code{\link{checkLog}},
#'    \code{\link{handleData}}, \code{\link{cureNegative}}, and \code{\link{remove_duplicate_genes}} as needed.
#'  \item **Build a \code{gene_networks} object** that includes the data matrix, sample metadata, pathway classifications, and the group and timepoints to be compared.
#'    \preformatted{
#'    net <- new("gene_network",
#'              gene_intensities = gi,                  #genes x samples matrix
#'              group_name       = 'Healthy',           #group name
#'              timepoints       = c(1,2),              #two time points
#'              metadata         = metadata df,
#'              output_dir       = 'path/to/results/',
#'              dataset          = 'GSEXXXXX',
#'              database         = 'KEGG',              #or GO
#'              pathways         = pathways,
#'              group_col        = "Group",             #group name in metadata
#'              subjects_ids     = "id",
#'              sample_ids       = "ID",
#'              time_col         = "Day",
#'              time_labels      = c('Day1', 'Day2))
#'    }
#'  \item **Run pipeline**:
#'    \preformatted{run_analyses(net)}
#'      This will call internally a) \code{\link{network_construction}} to build adjacency networks per time point and export edge lists
#'      b) \code{\link{compute_rewiring}} to estimate gene rewiring from one time point to another and export the rewiring scores per pathway (network)
#'      c) \code{\link{compute_statistics}} to get statistical summary e.g. median, mean, that results in one score per pathway as its overall rewiring.
#' }
#'
#' @section Output format:
#' \code{run_analyses} will export results in the directory \code{output_dir/dataset/database/<group>_<timepoint1>_<timepoint2>/}, with subfolders
#' \code{edge_lists/} (one edge list per pathway) and \code{dynet_score/} (rewiring scores per pathway), and a .RData file of the statistical summary.
#'
#' @keywords internal
"_PACKAGE"
