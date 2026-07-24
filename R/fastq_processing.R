# ==================================================
# FASTQ PROCESSING FUNCTIONS
# ==================================================
#
# Required packages:
# library(Biostrings)
# library(dplyr)
# library(data.table)
#
# ==================================================

# --------------------------------------------------
# Extract insert between flanking sequences
# --------------------------------------------------

extract_insert <- function(seq, left, right) {
  
  start <- regexpr(
    left,
    seq,
    fixed = TRUE
  )[1]
  
  if (start < 1)
    return(NA_character_)
  
  start <- start + nchar(left)
  
  end <- regexpr(
    right,
    substring(seq, start),
    fixed = TRUE
  )[1]
  
  if (end < 1)
    return(NA_character_)
  
  substring(
    seq,
    start,
    start + end - 2
  )
}

# --------------------------------------------------
# Try forward and reverse-complement orientation
# --------------------------------------------------

extract_both <- function(seq, left_flank, right_flank) {
  
  seq <- gsub("[^ACGT]", "", seq)
  
  if (nchar(seq) < 10)
    return(NA_character_)
  
  x <- extract_insert(seq, left_flank, right_flank)
  
  if (!is.na(x))
    return(x)
  
  rc <- as.character(
    Biostrings::reverseComplement(
      Biostrings::DNAString(seq)
    )
  )
  
  extract_insert(rc, left_flank, right_flank)
}

# --------------------------------------------------
# Translate DNA to protein
# --------------------------------------------------

translate_seq <- function(seq) {
  
  n <- nchar(seq)
  
  if (n < 3)
    return(NA_character_)
  
  trimmed <- substr(
    seq,
    1,
    n - (n %% 3)
  )
  
  tryCatch(
    as.character(
      Biostrings::translate(
        Biostrings::DNAString(trimmed)
      )
    ),
    error = function(e) NA_character_
  )
  
}

# --------------------------------------------------
# Needleman-Wunsch alignment
# --------------------------------------------------

align_pair <- function(a, b) {
  
  p <- strsplit(a, "")[[1]]
  r <- strsplit(b, "")[[1]]
  
  n <- length(p)
  m <- length(r)
  
  score <- matrix(0, n + 1, m + 1)
  
  gap <- -1
  match <- 1
  mismatch <- -1
  
  for (i in 1:(n + 1))
    score[i, 1] <- (i - 1) * gap
  
  for (j in 1:(m + 1))
    score[1, j] <- (j - 1) * gap
  
  for (i in 2:(n + 1)) {
    for (j in 2:(m + 1)) {
      
      score[i, j] <- max(
        score[i - 1, j - 1] +
          ifelse(p[i - 1] == r[j - 1], match, mismatch),
        
        score[i - 1, j] + gap,
        
        score[i, j - 1] + gap
      )
    }
  }
  
  i <- n + 1
  j <- m + 1
  
  pa <- c()
  ra <- c()
  
  while (i > 1 || j > 1) {
    
    if (
      i > 1 &&
      j > 1 &&
      score[i, j] ==
      score[i - 1, j - 1] +
      ifelse(p[i - 1] == r[j - 1], match, mismatch)
    ) {
      
      pa <- c(p[i - 1], pa)
      ra <- c(r[j - 1], ra)
      
      i <- i - 1
      j <- j - 1
      
    } else if (
      i > 1 &&
      score[i, j] == score[i - 1, j] + gap
    ) {
      
      pa <- c(p[i - 1], pa)
      ra <- c("-", ra)
      
      i <- i - 1
      
    } else {
      
      pa <- c("-", pa)
      ra <- c(r[j - 1], ra)
      
      j <- j - 1
    }
  }
  
  list(
    prot = paste(pa, collapse = ""),
    ref  = paste(ra, collapse = "")
  )
}

# --------------------------------------------------
# Mutation annotation
# --------------------------------------------------

annotate <- function(prot, ref, count = NULL, freq = NULL) {
  
  if (nchar(prot) == nchar(ref) && !grepl("\\*", prot)) {
    
    p <- strsplit(prot, "")[[1]]
    r <- strsplit(ref, "")[[1]]
    
    idx <- which(p != r)
    
    if (length(idx)) {
      
      muts <- paste0(
        r[idx],
        idx,
        p[idx]
      )
      
      pos_df <- lapply(
        muts,
        function(x) {
          data.frame(
            variant = x,
            count = count,
            freq = freq
          )
        }
      )
      
    } else {
      
      muts <- character()
      pos_df <- list()
      
    }
    
    label <- paste(muts, collapse = ";")
    
    if (label == "")
      label <- "WT"
    
    return(
      list(
        label = label,
        pos = if (length(pos_df))
          dplyr::bind_rows(pos_df)
        else
          NULL
      )
    )
  }
  
  if (prot == ref)
    return(list(label = "WT"))
  
  if (grepl("\\*", prot))
    return(list(label = "STOP"))
  
  aln <- align_pair(prot, ref)
  
  p <- strsplit(aln$prot, "")[[1]]
  r <- strsplit(aln$ref, "")[[1]]
  
  ref_pos <- 0
  muts <- character()
  pos_df <- list()
  
  for (i in seq_along(p)) {
    
    if (r[i] != "-")
      ref_pos <- ref_pos + 1
    
    if (p[i] == "-" && r[i] != "-") {
      
      muts <- c(
        muts,
        paste0(r[i], ref_pos, "del")
      )
      
      pos_df[[length(pos_df) + 1]] <- data.frame(
        variant = paste0(r[i], ref_pos, "del"),
        count = count,
        freq = freq
      )
    }
    
    if (
      p[i] != r[i] &&
      p[i] != "-" &&
      r[i] != "-"
    ) {
      
      muts <- c(
        muts,
        paste0(r[i], ref_pos, p[i])
      )
      
      pos_df[[length(pos_df) + 1]] <- data.frame(
        variant = paste0(r[i], ref_pos, p[i]),
        count = count,
        freq = freq
      )
    }
  }
  
  list(
    label = paste(muts, collapse = ";"),
    pos = if (length(pos_df))
      dplyr::bind_rows(pos_df)
    else
      NULL
  )
}

# --------------------------------------------------
# Variant classification
# --------------------------------------------------

classify_variant <- function(prot, ref) {
  
  if (prot == ref)
    return("WT")
  
  if (grepl("\\*", prot))
    return("frameshift")
  
  if (nchar(prot) == nchar(ref))
    return("missense")
  
  if (nchar(prot) < nchar(ref))
    return("inframe_deletion")
  
  "frameshift"
}

# --------------------------------------------------
# Collapse identical mutation labels
# --------------------------------------------------

collapse_variants <- function(df) {
  
  as.data.table(df)[
    ,
    .(
      count = sum(count),
      freq = sum(freq)
    ),
    by = mutation
  ][
    order(
      -as.integer(mutation == "WT"),
      -freq
    )
  ]
  
}

# ==================================================
# STEP 2 PIPELINE
# FASTQ -> HAPLOTYPE TABLE
# ==================================================

process_haplotypes <- function(
    fastq_file,
    ref_protein,
    left_flank,
    right_flank,
    min_count = 10,
    min_freq = 0
) {
  t_total <- Sys.time()
  message("Reading FASTQ...")
  
  
  lines <- readLines(fastq_file)
  
  seqs_raw <- lines[
    seq(2, length(lines), 4)
  ]
  
  rm(lines)
  
  message(
    "Total reads: ",
    length(seqs_raw)
  )
  t <- Sys.time()
  
  message("Extracting inserts...")
  
  t_extract <- Sys.time()
  
  n_cores <- max(
    1,
    parallel::detectCores() - 1
  )
  
  seqs <- vapply(
    seqs_raw,
    extract_both,
    character(1),
    left_flank = left_flank,
    right_flank = right_flank
  )
  message(
    "Insert extraction took ",
    round(
      as.numeric(
        difftime(
          Sys.time(),
          t_extract,
          units = "secs"
        )
      ),
      2
    ),
    " sec"
  )
  
  seqs <- unlist(
    seqs,
    use.names = FALSE
  )
  
  seqs <- seqs[!is.na(seqs)]
  
  message(
    "Valid inserts: ",
    length(seqs)
  )
  
  message("Counting haplotypes...")
  
  haplo_df <- data.table::data.table(
    dna = seqs
  )[
    ,
    .(count = .N),
    by = dna
  ]
  
  haplo_df[
    ,
    freq := count / sum(count)
  ]
  
  haplo_df_raw <- data.table::copy(
    haplo_df
  )
  
  haplo_df <- haplo_df[
    count >= min_count &
      freq >= min_freq
  ]
  t <- Sys.time()
  message(
    "Unique haplotypes (raw): ",
    nrow(haplo_df_raw)
  )
  
  message(
    "Unique haplotypes (filtered): ",
    nrow(haplo_df)
  )
  message(
    "Mean insert length: ",
    round(
      mean(
        nchar(haplo_df_raw$dna)
      ),
      1
    )
  )
  
  message("Translating proteins...")
  
  # Translate filtered table only
  
  haplo_df[
    ,
    protein := vapply(
      dna,
      translate_seq,
      character(1)
    )
  ]
  
  message(
    "Translation took ",
    round(
      as.numeric(
        difftime(
          Sys.time(),
          t,
          units = "secs"
        )
      ),
      2
    ),
    " sec"
  )
  
  qc <- list(
    total_reads = length(seqs_raw),
    valid_inserts = length(seqs),
    unique_haplotypes_raw = nrow(haplo_df_raw),
    unique_haplotypes_filtered = nrow(haplo_df)
  )
  message(
    "Total runtime: ",
    round(
      as.numeric(
        difftime(
          Sys.time(),
          t_total,
          units = "secs"
        )
      ),
      2
    ),
    " sec"
  )
  list(
    haplo_df = haplo_df,
    haplo_df_raw = haplo_df_raw,
    qc = qc
  )
  
}

# ==================================================
# HAPLOTYPE TABLE -> MUTATION COUNTS
# ==================================================
annotate_variants <- function(
    haplo_df,
    haplo_df_raw,
    ref_protein
) {
  
  annotate_table <- function(df, ref_protein) {
    
    ann <- mapply(
      annotate,
      df$protein,
      ref_protein,
      df$count,
      df$freq,
      SIMPLIFY = FALSE
    )
    
    df$mutation <- vapply(
      ann,
      `[[`,
      character(1),
      "label"
    )
    
    df$variant_class <- mapply(
      classify_variant,
      df$protein,
      ref_protein
    )
    
    df$n_mut <- lengths(
      strsplit(df$mutation, ";")
    )
    
    df
  }
  
  # Annotate filtered table only
  
  haplo_df <- annotate_table(
    haplo_df,
    ref_protein
  )
  
  single_df <- haplo_df[
    (haplo_df$n_mut == 1 |
       haplo_df$mutation == "WT") &
      haplo_df$mutation != "STOP",
  ]
  
  single_collapsed <- collapse_variants(single_df)
  
  missense_df <- collapse_variants(
    haplo_df[haplo_df$variant_class == "missense", ]
  )
  
  indel_df <- collapse_variants(
    haplo_df[haplo_df$variant_class == "inframe_deletion", ]
  )
  
  frameshift_df <- collapse_variants(
    haplo_df[haplo_df$variant_class == "frameshift", ]
  )
  
  list(
    haplo_df = haplo_df,
    haplo_df_raw = haplo_df_raw,
    single_mutants = single_collapsed,
    missense = missense_df,
    indels = indel_df,
    frameshifts = frameshift_df
  )
}
# ==================================================
# COMBINE INTO SINGLE WRAPPER FUNCTION
# ==================================================
run_mutation_analysis <- function(
    fastq_file,
    ref_protein,
    left_flank,
    right_flank,
    min_count,
    min_freq
) {
  
  haplo_res <- process_haplotypes(
    fastq_file = fastq_file,
    ref_protein = ref_protein,
    left_flank = left_flank,
    right_flank = right_flank,
    min_count = min_count,
    min_freq = min_freq
  )
  
  variant_res <- annotate_variants(
    haplo_df = haplo_res$haplo_df,
    haplo_df_raw = haplo_res$haplo_df_raw,
    ref_protein = ref_protein
  )
  
  variant_res$qc <- haplo_res$qc
  
  variant_res
}

