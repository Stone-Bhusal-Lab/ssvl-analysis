# ============================================================================
# Package Management
# ============================================================================

# ============================================================================
# PACKAGE DEFINITIONS
# ============================================================================

cran_packages <- c(
  "shiny",
  "dplyr",
  "stringr",
  "ggplot2",
  "plotly",
  "ggrepel",
  "circlize",
  "grid",
  "DT",
  "r3dmol",
  "htmlwidgets",
  "uuid",
  "tibble",
  "purrr",
  "tidyr",
  "readr",
  "data.table",
  "parallel"
)

bioc_packages <- c(
  "Biostrings",
  "ShortRead",
  "ComplexHeatmap"
)

required_packages <- c(
  cran_packages,
  bioc_packages
)

# ============================================================================
# CHECK REQUIRED PACKAGES
# ============================================================================

check_required_packages <- function() {
  
  installed <- rownames(
    installed.packages()
  )
  
  missing <- setdiff(
    required_packages,
    installed
  )
  
  if (length(missing) > 0) {
    
    stop(
      paste0(
        "\n\nMissing required packages:\n\n",
        paste(
          missing,
          collapse = "\n"
        ),
        "\n\nRun:\n",
        "source('package_setup.R')\n",
        "to install missing dependencies.\n"
      ),
      call. = FALSE
    )
    
  }
  
  message(
    "All required packages installed."
  )
  
}
# ============================================================================
# Installer function
# ============================================================================

install_missing_packages <- function() {
  
  installed <- rownames(installed.packages())
  
  cran_missing <- setdiff(
    cran_packages,
    installed
  )
  
  if (length(cran_missing) > 0) {
    
    install.packages(cran_missing)
    
  }
  
  if (!requireNamespace(
    "BiocManager",
    quietly = TRUE
  )) {
    
    install.packages("BiocManager")
    
  }
  
  bioc_missing <- setdiff(
    bioc_packages,
    rownames(installed.packages())
  )
  
  if (length(bioc_missing) > 0) {
    
    BiocManager::install(
      bioc_missing,
      ask = FALSE,
      update = FALSE
    )
    
  }
  
}
# ============================================================================
# Loader function
# ============================================================================

load_required_packages <- function() {
  
  invisible(
    lapply(
      c(cran_packages, bioc_packages),
      library,
      character.only = TRUE
    )
  )
  
}