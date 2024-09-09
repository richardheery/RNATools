#' Extract forward and reverse reads from BAM files and save to FASTQ files.
#' 
#' Took 4.5 minutes with a BAM of size 2.8 GB.
#'
#' @param path_to_samtools Path to samtools directory
#' @param bam_files A vector of paths to BAM files
#' @param fastq_file_prefixes Prefixes to use for output FASTQ files. 
#' "_f1.fq" is added to the prefixes for forward reads and "_r2.fq" is added for reverse reads. Default is the basenames of bam_files. 
#' @param output_directory Path to output directory to save files FASTQ files. Directory is created if it doesn't exist.  Default is current working directory. 
#' @param parallel_cores The number of files to process at a time in parallel. Default is just 1
#' @param compression One of ".gz" to compress files with gzip, ".bgzf" to compress files with bgzip or NULL to not compress files. Default is gzip compression.
#' @export
samtools_fastq = function(path_to_samtools, bam_files, fastq_file_prefixes = NULL, output_directory = ".", parallel_cores = 1, compression = ".gz"){
  
  # Check that one of the allowed values is given for compression
  match.arg(arg = compression, choices = c(".gz", ".bgzf", NULL), several.ok = F)
  
  # Get the canonical path to samtools
  path_to_samtools = normalizePath(path_to_samtools)
  
  # Check that samtools works
  if(system2(path_to_samtools, c("--version"), stdout = F, stderr = F) != 0){stop("samtools cannot be called from given path")}
  
  # Create fastq_file_prefixes if they are not provided
  if(is.null(fastq_file_prefixes)){
    fastq_file_prefixes = basename(tools::file_path_sans_ext(bam_files))
  }
  
  # Check that fastq_file_prefixes has the same length as bam_files
  if(length(bam_files) != length(fastq_file_prefixes)){
    stop("fastq_file_prefixes must have the same length as bam_files")
  }
  
  # Create output_directory if it doesn't exist
  if(!dir.exists(output_directory)){dir.create(output_directory)}
  
  # Name the forward and reverse fastq files. 
  forward_fastqs = paste0(output_directory, "/", fastq_file_prefixes, "_f1.fq", compression)
  reverse_fastqs = paste0(output_directory, "/", fastq_file_prefixes, "_r2.fq", compression)
  
  # Create cluster if parallel_cores greater than 1
  if(parallel_cores > 1){
    cl = parallel::makeCluster(parallel_cores)
    doParallel::registerDoParallel(cl, parallel_cores)
    `%dopar%` = foreach::`%dopar%`
    } else {`%dopar%` = foreach::`%do%`}
  
  foreach::foreach(bam = seq_along(bam_files)) %dopar% {
    system2(command = path_to_samtools, 
      args = c("fastq", "-1", forward_fastqs[bam], "-2", reverse_fastqs[bam],  
        "-0 /dev/null", "-s /dev/null", "-n", bam_files[bam]))
  }
}

#' Sort BAM files
#' 
#' @param path_to_samtools Path to samtools directory
#' @param bam_files Paths to BAM files
#' @param sorted_bam_files Paths to save sorted BAM files. Note that sortBam adds .bam to end of file even if it is already present. 
#' Default is to add "sorted" to end of names of bam_files
#' @param parallel_cores Number of cores to use. Default is 1.
#' @export
sort_bams = function(path_to_samtools, bam_files, sorted_bam_files = NULL, parallel_cores = 1){
  
  # Get the canonical path to samtools
  path_to_samtools = normalizePath(path_to_samtools)
  
  # Check that samtools works
  if(system2(path_to_samtools, c("--version"), stdout = F, stderr = F) != 0){stop("samtools cannot be called from given path")}
  
  # If sorted_bam_files not provided, add "sorted" to end of names of bam_files
  if(is.null(sorted_bam_files)){
    sorted_bam_files = gsub(".bam", "_sorted", bam_files)
  }
  
  # Check that bam_files and sorted_bam_files have the same length and that there are no duplications
  if(length(bam_files) != length(sorted_bam_files)){stop("bam_files and sorted_bam_files must have the same length")}
  if(anyDuplicated(bam_files) | anyDuplicated(sorted_bam_files)){
    stop("There cannot be any duplicated names in bam_files or sorted_bam_files")}
  
  `%dopar%` = foreach::`%do%`
  foreach::foreach(bam = seq_along(bam_files)) %dopar% {
    system2(command = path_to_samtools, 
      args = c("sort", 
        bam_files[bam],
        "-o", sorted_bam_files[bam],
        "-@", parallel_cores
      )
    )
  }
}


#' Index BAM files
#' 
#' @param path_to_samtools Path to samtools directory
#' @param bam_files Paths to BAM files
#' @param parallel_cores Number of cores to use. Default is 1.
#' @export
index_bams = function(path_to_samtools, bam_files, parallel_cores = 1){
  
  # Get the canonical path to samtools
  path_to_samtools = normalizePath(path_to_samtools)
  
  # Check that samtools works
  if(system2(path_to_samtools, c("--version"), stdout = F, stderr = F) != 0){stop("samtools cannot be called from given path")}
  
  `%dopar%` = foreach::`%do%`
  foreach::foreach(bam = seq_along(bam_files)) %dopar% {
    system2(command = path_to_samtools, 
      args = c("index", 
        bam_files[bam],
        "-@", parallel_cores
      )
    )
  }
}


#' Calculate coverage of BAM files using bamCoverage from deepTools
#' 
#' @param path_to_bamCoverage Path to bamCoverage
#' @param bam_files Paths to BAM files
#' @param output_bigwigs Paths to save output bigWigs. Default is to replace ".bam" with ".bw". 
#' @param parallel_cores Number of cores to use. Default is 1.
#' @export
bam_coverage = function(path_to_bamCoverage, bam_files, output_bigwigs = NULL, binSize = 50, effectiveGenomeSize, parallel_cores = 1){
  
  # If output_bigwigs not provided, change ".bam" to ".bw"
  if(is.null(output_bigwigs)){
    output_bigwigs = gsub(".bam", ".bw", bam_files)
  }
  
  # Get the canonical path to samtools
  path_to_bamCoverage = normalizePath(path_to_bamCoverage)
  
  # Check that bam_files and output_bigwigs have the same length and that there are no duplications
  if(length(bam_files) != length(output_bigwigs)){stop("bam_files and output_bigwigs must have the same length")}
  if(anyDuplicated(bam_files) | anyDuplicated(output_bigwigs)){
    stop("There cannot be any duplicated names in bam_files or output_bigwigs")}
  
  # Check that there is no overlap between bam_files and output_bigwigs
  if(length(intersect(bam_files, output_bigwigs))){stop("Some output_bigwigs overlap bam_files")}
  
  `%do%` = foreach::`%do%`
  foreach::foreach(bam = seq_along(bam_files)) %do% {
    print(paste0("Starting BAM file ", bam))
    system2(command = path_to_bamCoverage, 
      args = c(
        "-b", bam_files[bam], 
        "-o", output_bigwigs[bam], 
        "--effectiveGenomeSize", --effectiveGenomeSize,
        "-p", parallel_cores,
        "--binSize", binSize
        )
    )
  }
}



