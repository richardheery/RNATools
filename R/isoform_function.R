#' Quantify relative isoform abundance
#'
#' @param isoform_expression_table A table with the expression values of gene isoforms. 
#' Row names should be the names of the isoforms.
#' @param isoform_grl A GRangesList object where the name of each component GRanges 
#' is the name of a gene and the name of the individual regions are the isoform names. 
#' @param ties.method A character string specifying how ties are treated when ranking isoform expression. 
#' Same as the options for base::rank. Default behavior is to use the mean of tied ranks. 
#' @return A RangedSummarizedExperiment with the relative abundance and expression rank of isoforms.
#' @export
relative_isoform_abundance = function(isoform_expression_table, isoform_grl, 
  ties.method = c("average", "first", "last", "random", "max", "min")){
  
  # Check that one of the permitted values provided for ties.method
  ties.method = match.arg(ties.method)
  
  # Create a data.frame matching genes to isoforms
  gene_to_isoform = data.frame(
    gene = rep(names(isoform_grl), times = lengths(isoform_grl)),
    isoform = names(unlist(unname(isoform_grl)))
  )
  
  # Create a flat GRanges from isoform_grl and add gene name as a metadata column
  isoform_gr_flat = unlist(unname(isoform_grl))
  mcols(isoform_gr_flat) = NULL
  isoform_gr_flat$gene = gene_to_isoform$gene[match(names(isoform_gr_flat), gene_to_isoform$isoform)]
  
  # Convert row.names of isoform_expression_table to a column called isoform
  isoform_expression_table = tibble::rownames_to_column(isoform_expression_table, "isoform")
  
  # Convert isoform_expression_table to long format
  isoform_expression_table_long = tidyr::pivot_longer(isoform_expression_table, 
    -isoform, names_to = "sample", values_to = "expression")
  
  # Add a column with gene name
  isoform_expression_table_long$gene = gene_to_isoform$gene[
    match(isoform_expression_table_long$isoform, gene_to_isoform$isoform)]
  
  # Group isoform_expression_table_long by gene and sample
  isoform_expression_table_long = dplyr::group_by(isoform_expression_table_long, gene, sample)
  
  # Add the relative abundance and expression ranking of each isoform
  isoform_expression_table_long = dplyr::mutate(isoform_expression_table_long, 
    relative_abundance = expression/sum(expression))
  isoform_expression_table_long = dplyr::mutate(isoform_expression_table_long, 
    isoform_rank = rank(-relative_abundance, ties.method = ties.method))
  
  # If relative_abundance is NaN, set isoform_rank to NaN
  isoform_expression_table_long$isoform_rank[is.na(isoform_expression_table_long$relative_abundance)] = NaN
  
  # Ungroup isoform_expression_table_long
  isoform_expression_table_long = dplyr::ungroup(isoform_expression_table_long)
  
  # Create a table with relative isoform abundance in wide format
  isoform_abundance_table = tidyr::pivot_wider(
    dplyr::select(isoform_expression_table_long, isoform, sample, relative_abundance), 
    names_from = "sample", values_from = "relative_abundance")
  
  # Set isoform to row.names and put rows in order of isoforms in isoform_gr_flat
  isoform_abundance_table = tibble::column_to_rownames(isoform_abundance_table, "isoform")
  isoform_abundance_table = isoform_abundance_table[names(isoform_gr_flat), ]
  
  # Create a table with isoform rank in wide format
  isoform_rank_table = tidyr::pivot_wider(
    dplyr::select(isoform_expression_table_long, isoform, sample, isoform_rank), 
    names_from = "sample", values_from = "isoform_rank")
  
  # Set isoform to row.names and put rows in order of isoforms in isoform_gr_flat
  isoform_rank_table = tibble::column_to_rownames(isoform_rank_table, "isoform")
  isoform_rank_table = isoform_rank_table[names(isoform_gr_flat), ]
  
  # Create a RangedSummarizedExperiment with the isoform abundance and rank tables
  isoform_se = SummarizedExperiment::SummarizedExperiment(
    assays = list("relative_abundance" = isoform_abundance_table, "isoform_rank" = isoform_rank_table),
    rowRanges = isoform_gr_flat)
  
}

#' Combine the expression values of isoforms to get overall expression of their associated genes
#' 
#' @param isoform_expression_table A table where rows are isoforms and columns are samples. Row names should be the names of isoforms. 
#' @param gene_to_isoform_list A list of vectors where the name of each list entry is a gene name and its elements are the names of isoforms.
#' Can alternatively be a GRangeList where the name of each list element is a gene and the names of the individual ranges are isoforms.
#' @return A data.frame with the sum of isoform expression values for genes where rows are genes and columns are samples
#' @export
sum_isoform_values <- function(isoform_expression_table, gene_to_isoform_list){
  
  # Check that inputs have the correct data type
  stopifnot(is(isoform_expression_table, "data.frame") | is(isoform_expression_table, "matrix"), 
    is(gene_to_isoform_list, "list"), all(sapply(gene_to_isoform_list, function(x) is(x, "character"))))
  
  # If gene_to_isoform_list is a GRangesList, extract a list of vectors matching genes to isoforms
  if(is(gene_to_isoform_list, "GRangesList")){
    gene_to_isoform_list <- lapply(gene_to_isoform_list, names)
  }
  
  # Get the sum of the expression values for all isoforms associated with each gene in each sample
  `%do%` <- foreach::`%do%`
  results_list <- foreach::foreach(gene_isoforms = gene_to_isoform_list) %do% {
    
    # Sum the isoform values for each isoform associated with a gene
    gene_results <- colSums(isoform_expression_table[gene_isoforms, ], na.rm = TRUE)
    
    # Samples where all values for gene_isoforms are NA are given a value of NA for the gene. 
    # This is done as colSums returns a value of 0 if all values in a column are NA and na.rm <- TRUE.  
    gene_results[apply(isoform_expression_table[gene_isoforms, ], 2, function(x) all(is.na(x)))] <- NA
    gene_results
  }
  
  # Set names for results_list
  names(results_list) <- names(gene_to_isoform_list)
  
  # Turn results_list into a data.frame
  results_table <- data.frame(dplyr::bind_rows(results_list, .id = "gene_name"))
  
  # Set gene names as row.names and return
  return(tibble::column_to_rownames(results_table, "gene_name"))
}