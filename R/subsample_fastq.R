head_fastq <- function(
    input_fastq,
    output_fastq,
    n_reads = 5000
) {
  
  n_lines <- n_reads * 4
  
  con_in <- if (
    grepl("\\.gz$", input_fastq)
  ) {
    gzfile(input_fastq, "rt")
  } else {
    file(input_fastq, "r")
  }
  
  on.exit(close(con_in))
  
  lines <- readLines(
    con_in,
    n = n_lines
  )
  
  writeLines(
    lines,
    output_fastq
  )
  
}

