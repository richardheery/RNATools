#' Run salmon on a group of FASTQ files
#'
#' @param path_to_salmon Path to Salmon executable
#' @param salmon_index Path to a salmon index 
#' @param forward_fastqs A vector with the paths to forward FASTQ files. Each file should be the match of the corresponding file in reverse_fastqs.
#' @param reverse_fastqs A vector with the paths to reverse FASTQ files. Each file should be the match of the corresponding file in forward_fastqs.
#' @param sample_names A vector with the names of samples corresponding for each pair of samples from forward_fastqs and reverse_fastqs
#' @param output_directory The name of the directory to save results in. Will be created if it doesn't already exist. 
#' @param merged_output_prefix Prefix to use for names of merged output files for counts and TPM which take the form 
#' merged_output_prefix_counts_merged.tsv.gz and merged_output_prefix_tpm_merged.tsv.gz. Default prefix is "salmon_transcript" i.e. default output merged output files are 
#' salmon_transcript_counts_merged.tsv.gz and salmon_transcript_tpm_merged.tsv.gz. 
#' @param merged_files The names of the files in output_directory to save the count and TPM values for all samples to. Defaults are salmon_transcript_count_merged.tsv and salmon_transcript_tpm_merged.tsv
#' @param messages_file Name of file to save Salmon run messages. If no file name given, information is printed to stdout.
#' @param n_cores The number of cores to use. Default is 1.
#' @export
salmon_quantify = function(path_to_salmon, salmon_index, forward_fastqs, reverse_fastqs, 
  sample_names, output_directory, merged_output_prefix = "salmon_transcript", messages_file = "", n_cores = 1){
  
  # Check if salmon can be executed from the given path
  if(suppressWarnings(system(paste(path_to_salmon, "-h"), ignore.stdout = T, ignore.stderr = F)) != 0){stop("Salmon could can not be executed from given path")}
  
  # Get the canonical path to Salmon
  path_to_salmon = normalizePath(path_to_salmon)
  
  # Check that path to salmon_index exits
  if(!dir.exists(salmon_index)){stop("salmon_index could not be found")}
  
  # Check that length of forward_fastqs and reverse_fastqs
  if(length(forward_fastqs) != length(reverse_fastqs)){stop("forward_fastqs and reverse_fastqs have different lengths")}
  if(length(forward_fastqs) != length(sample_names)){stop("sample_names does not have same length as forward_fastqs and reverse_fastqs")}
  
  # Check if messages_file already exists
  if(messages_file != ""){if(file.exists(messages_file)){stop(paste(messages_file, "already exists"))}}
  
  # Create output_directory if it doesn't exist
  if(!dir.exists(output_directory)){dir.create(output_directory)}
  
  # Create output subdirectories for each pair of FASTQ samples by pasting sample_names to output_directory
  sample_directories = paste(output_directory, sample_names, sep = "/")
  
  # Check that sample directories and merged_files do not already exist
  for(directory in sample_directories) {
    if(file.exists(directory)){stop(paste(directory, "already exists"))}
  }
  
  # Create names for merged output files and check that they do not 
  merged_counts_file = paste(merged_output_prefix, "counts_merged.tsv.gz", sep = "_")
  merged_tpm_file = paste(merged_output_prefix, "tpm_merged.tsv.gz", sep = "_")
  
  if(file.exists(paste(output_directory, merged_counts_file, sep = "/"))){
    stop(paste(paste(output_directory, merged_counts_file, sep = "/"), "already exists"))}
  
  if(file.exists(paste(output_directory, merged_tpm_file, sep = "/"))){
    stop(paste(paste(output_directory, merged_tpm_file, sep = "/"), "already exists"))}
  
  # Loop through each pair of FASTQs and run salmon
  for(pair in seq_along(forward_fastqs)){
    print(paste("Starting to process FASTQ pair", pair))
    system2(command = path_to_salmon,
      args = paste("quant -i", salmon_index, "-l A", "-1", forward_fastqs[pair], "-2", reverse_fastqs[pair], "-p", n_cores, "-o", sample_directories[pair]), stderr = messages_file)
  }
  
  # Combine TPM from output into a single file
  system2(command = path_to_salmon, paste("quantmerge --quants", paste(sample_directories, collapse = " "), "--column numreads", "-o",  paste(output_directory, merged_counts_file, sep = "/")), stdout = NULL, stderr = NULL)
  system2(command = path_to_salmon, paste("quantmerge --quants", paste(sample_directories, collapse = " "), "-column TPM", "-o",  paste(output_directory, merged_tpm_file, sep = "/")), stdout = NULL, stderr = NULL)
}