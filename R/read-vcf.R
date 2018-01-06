#' Read a VCF table
#'
#' @param file `[character(1)]` Path to VCF file. Supported file extensions:
#'   `.vcf`, `.vcf.gz`, `.vcf.bz2`. Compressed files are unzipped for the
#'   duration of the R session via [R.utils::gunzip], or [R.utils::bunzip2].
#'
#' @param col_types `[col_spec|character(1)|NULL]` Either NULL, or a valid
#'   column specification (for more details, see [readr::read_tsv]). Defaults
#'   to the following expected VCF columns:
#'
#'     - CHROM `[character]`
#'     - POS `[integer]`
#'     - ID `[character]`
#'     - REF `[character]`
#'     - ALT `[character]`
#'     - QUAL `[double]`
#'     - FILTER `[character]`
#'     - INFO `[character]`
#'
#' @param callback `[function(x, pos)|NULL]` A function with arguments `x` and
#'   `pos` that takes a dataframe as `x` and returns a dataframe (`pos` is used
#'   internally). This is usefull for filtering or summarizing VCF files that
#'   are otherwise too big to read into memory. See [readr::DataFrameCallback]
#'   and [readr::read_tsv_chunked] for examples. Note that when using a
#'   callback, the default value of `chunk_size` for [readr::read_tsv_chunked]
#'   has been increased to `1e6`. Defaults to `NULL`.
#'
#' @param ... Arguments passed to [readr::read_tsv], or
#'   [readr::read_tsv_chunked] if `callback` is not `NULL`.
#'
#' @return An object with class `vcf` (inherits `data.frame` class). VCF header
#'   information is stored as an attribute accessible via `attr(x, 'vcf')`.
#'
#' @references
#'   - https://faculty.washington.edu/browning/beagle/intro-to-vcf.html
#'
#' @importFrom tools file_ext
#' @importFrom R.utils gunzip bunzip2
#' @importFrom stringr str_split str_replace
#' @importFrom magrittr %>%
#' @import readr
#' @export
#' @md

read_vcf <- function(file,
                     col_types = readr::cols_only(
                       CHROM  = readr::col_character(),
                       POS    = readr::col_integer(),
                       ID     = readr::col_character(),
                       REF    = readr::col_character(),
                       ALT    = readr::col_character(),
                       QUAL   = readr::col_double(),
                       FILTER = readr::col_character(),
                       INFO   = readr::col_character()
                     ),
                     callback = NULL,
                     ...) {

  # Temporarily extract gz or bz2 compressed files
  file <-
    switch(
      tools::file_ext(file),
      vcf = file,
      gz  = R.utils::gunzip(file,  temporary = TRUE, skip = TRUE, remove = FALSE),
      bz2 = R.utils::bunzip2(file, temporary = TRUE, skip = TRUE, remove = FALSE)
    )

  # Locate header
  nline  = 0L
  header = integer(0)
  while (!length(header)) {
    lines  <- readr::read_lines(file, skip = nline, n_max = 300L)
    header <- grep('^#CHROM', lines) + nline
    nline  <- nline + 300L
  }

  # Extract column names
  col_names <-
    readr::read_lines(file, skip = header - 1L, n_max = 1L) %>%
    stringr::str_split('\t') %>%
    unlist() %>%
    stringr::str_replace('#CHROM', 'CHROM')

  # Read VCF table
  if (is.null(callback)) {
    vcf <-
      readr::read_tsv(
        file,
        skip      = header,
        na        = '.',
        col_names = col_names,
        col_types = col_types,
        ...
      )
  } else {

    if (missing(chunk_size)) chunk_size <- 1e6

    vcf <-
      readr::read_tsv_chunked(
        file,
        skip       = header,
        na         = '.',
        col_names  = col_names,
        col_types  = col_types,
        chunk_size = chunk_size,
        callback   = readr::DataFrameCallback$new(callback),
        ...
      )
  }

  # Append metadata as an attribute
  attr(vcf, 'vcf') <- readr::read_lines(file, n_max = header - 1L)
  class(vcf) <- c('vcf', class(vcf))
  return(vcf)
}
