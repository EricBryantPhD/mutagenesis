#' Read and write VCF tables
#'
#' @param file `[character(1)]` Path to VCF file. Supported file extensions:
#'   `.vcf`, `.vcf.gz`, `.vcf.bz2`. Compressed files are unzipped for the
#'   duration of the R session via [R.utils::gunzip], or [R.utils::bunzip2].
#'   `write_vcf` only supports `.vcf`.
#'
#' @param col_types `[col_spec|character(1)|NULL]` Either NULL, or a valid
#'   column specification (for more details, see [readr::read_tsv]). Defaults
#'   to only read the following expected VCF columns:
#'
#'     - CHROM [character]
#'     - POS [integer]
#'     - ID [character]
#'     - REF [character]
#'     - ALT [character]
#'     - QUAL [double]
#'     - FILTER [character]
#'     - INFO [character]
#'
#' @param callback `[function(x, pos)|NULL]` A function with arguments `x` and
#'   `pos` that takes a dataframe as `x` and returns a dataframe (`pos` is used
#'   internally). This is usefull for filtering or summarizing VCF files that
#'   are otherwise too big to read into memory. See [readr::DataFrameCallback]
#'   and [readr::read_tsv_chunked] for examples. Defaults to `NULL`.
#'
#' @param chunk_size `[integer(1)]` Number of lines to read for each chunk
#'   when using a callback. Defaults to `1e6`.
#'
#' @param ... Arguments passed to [readr::read_tsv], or
#'   [readr::read_tsv_chunked] if `callback` is not `NULL`.
#'
#' @param vcf `[vcf]` A VCF data frame. Simply a data frame with a `vcf` class
#'   and an attribute that contains VCF header information.
#'
#' @return An object with class `vcf` (inherits `data.frame` class). VCF header
#'   information is stored as an attribute accessible via `attr(x, 'vcf')`.
#'
#' @references
#'   - [Introduction to VCF](https://faculty.washington.edu/browning/beagle/intro-to-vcf.html)
#'
#' @examples
#'
#' @importFrom tools file_ext
#' @importFrom R.utils gunzip bunzip2
#' @importFrom stringr str_split str_replace
#' @importFrom magrittr %>%
#' @import readr
#' @export
#' @md

read_vcf <- function(file,
                     col_types = NULL,
                     callback = NULL,
                     chunk_size = 1e6,
                     ...) {

  if (is.null(col_types)) {
    col_types <- readr::cols_only(
      CHROM  = readr::col_character(),
      POS    = readr::col_integer(),
      ID     = readr::col_character(),
      REF    = readr::col_character(),
      ALT    = readr::col_character(),
      QUAL   = readr::col_double(),
      FILTER = readr::col_character(),
      INFO   = readr::col_character()
    )
  }

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
        na        = c('NA', '.'),
        col_names = col_names,
        col_types = col_types,
        ...
      )
  } else {
    vcf <-
      readr::read_tsv_chunked(
        file,
        skip       = header,
        na         =  c('NA', '.'),
        col_names  = col_names,
        col_types  = col_types,
        chunk_size = chunk_size,
        callback   = readr::DataFrameCallback$new(callback),
        ...
      )
  }

  # Convert missing ALT (i.e. '.') to empty string
  vcf <- mutate(vcf, ALT = dplyr::if_else(is.na(ALT), '', ALT))

  # Append metadata as an attribute
  attr(vcf, 'vcf') <- readr::read_lines(file, n_max = header - 1L)
  class(vcf) <- c('vcf', class(vcf))
  return(vcf)
}

#' @rdname read_vcf
#' @export
write_vcf <- function(vcf, file) {

  # Get header lines
  vcf_header <- attr(vcf, 'vcf')
  tbl_header <- paste0('#', paste(names(vcf), collapse = '\t'))

  # Convert missing ALT back to '.'
  vcf <- mutate(vcf, ALT = dplyr::if_else(ALT == '', '.', ALT))

  # Write to file
  readr::write_lines(vcf_header, file)
  readr::write_lines(tbl_header, file, append = TRUE)
  readr::write_tsv(vcf, file, na = '.', append = TRUE)
}
