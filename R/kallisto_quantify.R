#' Run kallisto on a group of FASTQ files
#'
#' @param path_to_kallisto Path to kallisto executable
#' @param kallisto_index Path to a kallisto index 
#' @param forward_fastqs A vector with the paths to forward FASTQ files. Each file should be the match of the corresponding file in reverse_fastqs.
#' @param reverse_fastqs A vector with the paths to reverse FASTQ files. Each file should be the match of the corresponding file in forward_fastqs.
#' @param sample_names A vector with the names of samples corresponding for each pair of samples from forward_fastqs and reverse_fastqs
#' @param output_directory The name of the directory to save results in. Is created if it doesn't already exist.
#' @param merged_output_prefix Prefix to use for names of merged output files for counts and TPM which take the form 
#' {merged_output_prefix}_counts_merged.tsv.gz and {merged_output_prefix}_tpm_merged.tsv.gz. Default prefix is "kallisto_transcript" i.e. default output merged output files are 
#' kallisto_transcript_counts_merged.tsv.gz and kallisto_transcript_tpm_merged.tsv.gz.
#' @param messages_file Name of file to save kallisto run messages. If no file name given, information is printed to stdout.
#' @param n_cores The number of cores to use. Default is 1.
#' @param number_bootstraps The number of bootstrap samples. Default is 100. 
#' @export
kallisto_quantify = function(path_to_kallisto, kallisto_index, forward_fastqs, reverse_fastqs, 
  sample_names, output_directory, merged_output_prefix = "kallisto_transcript", messages_file = "", n_cores = 1, number_bootstraps  = 100){
  
  # Check if kallisto can be executed from the given path
  if(suppressWarnings(system(paste(path_to_kallisto, "version"), ignore.stdout = T, ignore.stderr = F)) != 0){stop("kallisto could can not be executed from given path")}
  
  # Get the canonical path to kallisto
  path_to_kallisto = normalizePath(path_to_kallisto)
  
  # Check that path to kallisto_index exits
  if(!file.exists(kallisto_index)){stop("kallisto_index could not be found")}
  
  # Check that length of forward_fastqs, reverse_fastqs and sample names have the same lengths
  if(length(forward_fastqs) != length(reverse_fastqs)){stop("forward_fastqs and reverse_fastqs have different lengths")}
  if(length(forward_fastqs) != length(sample_names)){stop("sample_names does not have same length as forward_fastqs and reverse_fastqs")}
  
  # Check if messages_file already exists
  if(messages_file != ""){if(file.exists(messages_file)){stop(paste(messages_file, "already exists"))}}
  
  # Create output_directory if it doesn't exist
  if(!dir.exists(output_directory)){dir.create(output_directory)}
  
  # Create output subdirectories for each pair of FASTQ samples by pasting sample_names to output_directory
  sample_directories = paste(output_directory, sample_names, sep = "/")
  
  # Check that sample directories do not already exist
  for(directory in sample_directories) {
    if(file.exists(directory)){stop(paste(directory, "already exists"))}
  }
  
  # Create names for merged output files and check that they do not already exist
  merged_counts_file = paste(merged_output_prefix, "counts_merged.tsv.gz", sep = "_")
  merged_tpm_file = paste(merged_output_prefix, "tpm_merged.tsv.gz", sep = "_")
  
  if(file.exists(paste(output_directory, merged_counts_file, sep = "/"))){
    stop(paste(paste(output_directory, merged_counts_file, sep = "/"), "already exists"))}
  
  if(file.exists(paste(output_directory, merged_tpm_file, sep = "/"))){
    stop(paste(paste(output_directory, merged_tpm_file, sep = "/"), "already exists"))}
  
  # Loop through each pair of FASTQs and run kallisto
  for(pair in seq_along(forward_fastqs)){
    print(paste("Starting to process FASTQ pair", pair))
    system2(command = path_to_kallisto,
      args = paste("quant -i", kallisto_index, "-t", n_cores, "-b", number_bootstraps, "-o", sample_directories[pair], forward_fastqs[pair], reverse_fastqs[pair]),
      stderr = messages_file)
  }
  
  # Get paths to all abundance files
  abundance_files = paste0(sample_directories, "/abundance.tsv")
  
  # Create a data.frame with the counts calculated using kallisto for each sample
  kallisto_counts = data.frame(setNames(lapply(abundance_files, function(x) 
    data.table::fread(x, sep = "\t", header = T)$est_counts), basename(sample_directories)))
  
  # Create a data.frame with the TPM values calculated using kallisto for each sample
  kallisto_tpm = data.frame(setNames(lapply(abundance_files, function(x) 
    data.table::fread(x, sep = "\t", header = T)$tpm), basename(sample_directories)))
  
  # Get names of transcripts and add to output tables
  transcript_names = data.table::fread(abundance_files[1], sep = "\t", header = T)$target_id
  row.names(kallisto_counts) = transcript_names
  kallisto_counts = tibble::rownames_to_column(kallisto_counts, "transcript_id")
  row.names(kallisto_tpm) = transcript_names
  kallisto_tpm = tibble::rownames_to_column(kallisto_tpm, "transcript_id")
      
  # Write tables to output directory
  data.table::fwrite(kallisto_counts, paste(output_directory, merged_counts_file, sep = "/"), sep = "\t", row.names = F, quote = F)
  data.table::fwrite(kallisto_tpm, paste(output_directory, merged_tpm_file, sep = "/"), sep = "\t", row.names = F, quote = F)
  
}
