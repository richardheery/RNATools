#' Download SRA files with prefetch
#'
#' @param path_to_sratk Path to SRA toolkit bin directory
#' @param srr_accessions A vector of sequence read accessions
#' @param output_directory Path to output directory to save files. Default is current working directory. 
#' @param parallel_download The number of files to download at a time in parallel. Default is just 1
#' @param ngc_key Path to an ngc repository key. 
#' @export
sra_prefetch = function(path_to_sratk, srr_accessions, output_directory = ".", parallel_files = 1, ngc_key = NULL){
  
  # Get the canonical path to SRA toolkit
  path_to_sratk = normalizePath(path_to_sratk)
  path_to_prefetch = paste(path_to_sratk, "prefetch", sep = "/")
  
  # Check that prefetch works
  if(system2(path_to_prefetch, c("-V"), stdout = F, stderr = F) != 0){stop("prefetch cannot be called from given path")}
  
  if(!is.null(ngc_key)){
    ngc_flag = paste("--ngc", ngc_key)
  } else {
    ngc_flag = NULL
  }
  
  # Create cluster if parallel_files greater than 1
  cl = parallel::makeCluster(parallel_files)
  doParallel::registerDoParallel(cl, parallel_files)
  `%do%` = foreach::`%do%`
  `%dopar%` = foreach::`%dopar%`
  
  foreach::foreach(accession = srr_accessions) %dopar% {
    system2(command = path_to_prefetch, args = c(ngc_flag, paste("-O", output_directory), accession))
  }
}

#' Download SRA files with prefetch
#'
#' @param path_to_sratk Path to SRA toolkit bin directory
#' @param srr_directory_list A vector of sequence read accessions
#' @param output_directory Path to output directory to save files. Default is current working directory. 
#' @param parallel_download The number of files to download at a time in parallel. Default is just 1
#' @param remove_srr_directories A logical value indicating whether SRR directories should be removed after FASTQ files have been extracted. Default is FALSE. 
#' @export
sra_fastq_dump = function(path_to_sratk, srr_directory_list, output_directory = ".", parallel_files = 1, compress_fastq_files = T, remove_srr_directories = F){
  
  # Check that all SRR directories exist
  for(srr in srr_directory_list){
    if(!dir.exists(srr)){stop(paste(srr, "can't be found"))}
  }
  
  # Get the canonical path to SRA toolkit
  path_to_sratk = normalizePath(path_to_sratk)
  path_to_fastq_dump = paste(path_to_sratk, "fastq-dump", sep = "/")
  
  # Check that prefetch works
  if(system2(path_to_fastq_dump, c("-V"), stdout = F, stderr = F) != 0){stop("prefetch cannot be called from given path")}
  
  # Create cluster if parallel_files greater than 1
  cl = parallel::makeCluster(parallel_files)
  doParallel::registerDoParallel(cl, parallel_files)
  `%do%` = foreach::`%do%`
  `%dopar%` = foreach::`%dopar%`
  
  foreach::foreach(srr_directory = srr_directory_list) %dopar% {
    system2(command = path_to_fastq_dump, args = c("--split-3", paste("-O", output_directory), srr_directory))
    if(compress_fastq_files){
      for(fastq in list.files(path = output_directory, pattern = paste0(basename(srr_directory), ".*\\.fastq"), full.names = T)){
        R.utils::gzip(fastq)
      }
    }
    if(remove_srr_directories){
      unlink(srr_directory, recursive = T)
    }
  }
}