head_fastq <- function(
    input_fastq,
    output_fastq,
    n_reads = 5000
) {
  
  n_lines <- n_reads * 4
  
  con <- file(input_fastq, open = "r")
  
  on.exit(close(con))
  
  lines <- readLines(
    con,
    n = n_lines
  )
  
  writeLines(
    lines,
    output_fastq
  )
  
}