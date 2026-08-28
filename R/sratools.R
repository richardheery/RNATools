#' Download SRA files with SRA Toolkit prefetch
#'
#' @param srr_accessions A vector of sequence read accessions.
#' @param prefetch_path Optional path to SRA toolkit prefetch executable. Will attempt to search for prefetch on the system PATH if not provided.  
#' @param output_directory Path to output directory to save files. Default is current working directory. 
#' @param parallel_files The number of files to download at a time in parallel. Default is just 1.
#' @param ngc_key Path to an NGC repository key. 
#' @export
sra_prefetch = function(srr_accessions, prefetch_path = NULL, output_directory = ".", parallel_files = 1, ngc_key = NULL){
  
  # If prefetch_path is NULL, check if it is on the path
  if(is.null(prefetch_path)){
    prefetch_path = Sys.which("prefetch")
  }
  if(prefetch_path == ""){
    stop("Could not find prefetch on system PATH and prefetch_path is not provided")
  }
  
  # Check that prefetch works
  if(system2(prefetch_path, c("-V"), stdout = F, stderr = F) != 0){stop("prefetch cannot be called from ", prefetch_path)}
  
  # Define arguments to provide to prefetch
  prefetch_args = c("-O", output_directory)
  
  # If ngc_key provided, add to prefetch_args
  if(!is.null(ngc_key)){
    prefetch_args = c("--ngc", ngc_key, prefetch_args)
  }
  
  # Create cluster if parallel_files greater than 1
  if(parallel_files > 1){
    cl = parallel::makeCluster(parallel_files)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    doParallel::registerDoParallel(cl, parallel_files)
    `%dopar%` = foreach::`%dopar%`
  } else {
    `%dopar%` = foreach::`%do%`
  }
  
  # Prefetch files
  foreach::foreach(accession = srr_accessions) %dopar% {
    system2(command = prefetch_path, args = c(prefetch_args, accession))
  }
}

#' Extract prefetched SRA files with fastq-dump
#'
#' @param srr_directory_paths A vector of paths to SRR directories downloaded with prefetch.
#' @param fastq_dump_path Optional path to SRA toolkit fastq-dump executable. Will attempt to search for fastq-dump on the system PATH if not provided. 
#' @param output_directory Path to output directory to save files. Default is current working directory. 
#' @param parallel_files The number of files to extract at a time in parallel. Default is just 1.
#' @param compress_fastq_files Logical value indicating whether to compress FASTQ files with gzip after extracting. Default is TRUE.
#' @export
sra_fastq_dump = function(srr_directory_paths, fastq_dump_path = NULL, output_directory = ".", 
  parallel_files = 1, compress_fastq_files = T){
  
  # Check that all SRR directories exist
  for(srr in srr_directory_paths){
    if(!dir.exists(srr)){stop("Directory ", srr, " can't be found")}
  }
  
  # If fastq_dump_path is NULL, check if it is on the path
  if(is.null(fastq_dump_path)){
    fastq_dump_path = Sys.which("fastq-dump")
  }
  if(fastq_dump_path == ""){
    stop("Could not find fastq-dump on system PATH and fastq_dump_path is not provided")
  }
  
  # Check that fastq-dump works
  if(system2(fastq_dump_path, c("-V"), stdout = F, stderr = F) != 0){stop("fastq-dump cannot be called from ", fastq_dump_path)}
  
  # Create cluster if parallel_files greater than 1
  if(parallel_files > 1){
    cl = parallel::makeCluster(parallel_files)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    doParallel::registerDoParallel(cl, parallel_files)
    `%dopar%` = foreach::`%dopar%`
  } else {
    `%dopar%` = foreach::`%do%`
  }
  
  # Define arguments for fastq-dump
  fastq_dump_args = c("--split-3", "-O", output_directory)
  
  # Extract FASTQ files and compress with gzip if specified
  foreach::foreach(srr_directory = srr_directory_paths) %dopar% {
    system2(command = fastq_dump_path, args = c(fastq_dump_args, srr_directory))
    if(compress_fastq_files){
      for(fastq in list.files(path = output_directory, pattern = paste0("^", basename(srr_directory), ".*\\.fastq$"), full.names = T)){
        R.utils::gzip(fastq)
      }
    }
  }
}
