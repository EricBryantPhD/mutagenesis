#' Reverse and/or complement strings
#'
#' Given a character vector, returns the reverse, complement, or reverse &
#' complement of each string as a character vector.
#'
#' @param x `[character]`
#'
#'   Strings to reverse and/or complement.
#'
#' @param from_to `[named:character]`
#'
#'   Names and values correspond to letters and their complement. Defaults
#'   to `DNA`.
#'
#' @return
#'   Returns a character vector with the same length as `x`.
#'
#' @format
#'   **DNA** - A named `character` vector. Names and values
#'   correspond to complementary DNA bases.
#'
#'   **RNA** - A named `character` vector. Names and values
#'   correspond to complementary RNA bases.
#'
#' @examples
#' dna <- c('ATGATGC', 'AAAGGG')
#' reverse(dna)
#' complement(dna)
#' reverse_complement(dna)
#'
#' rna <- c('AUGAUGC', 'AAAGGG')
#' reverse(rna)
#' complement(rna, RNA)
#' reverse_complement(rna, RNA)
#'
#' @importFrom stringr str_c
#' @importFrom stringi stri_reverse
#' @export
#' @md

reverse_complement <- function(x, from_to = DNA) {
  reverse(complement(x, from_to))
}

#' @rdname reverse_complement
#' @export

reverse <- function(x) {
  stringi::stri_reverse(x)
}

#' @rdname reverse_complement
#' @export

complement <- function(x, from_to = DNA) {
  from <- stringr::str_c(names(from_to), collapse = '')
  to   <- stringr::str_c(from_to, collapse = '')
  chartr(from, to, x)
}

#' @rdname reverse_complement
#' @export

DNA <-
  c(
    'A' = 'T',
    'T' = 'A',
    'G' = 'C',
    'C' = 'G',
    'Y' = 'R',
    'R' = 'Y',
    'S' = 'S',
    'W' = 'W',
    'K' = 'M',
    'M' = 'K',
    'B' = 'V',
    'D' = 'H',
    'H' = 'D',
    'V' = 'B',
    'N' = 'N',
    'a' = 't',
    't' = 'a',
    'g' = 'c',
    'c' = 'g',
    'y' = 'r',
    'r' = 'y',
    's' = 's',
    'w' = 'w',
    'k' = 'm',
    'm' = 'k',
    'b' = 'v',
    'd' = 'h',
    'h' = 'd',
    'v' = 'b',
    'n' = 'n'
  )

#' @rdname reverse_complement
#' @export

RNA <-
  c(
    'A' = 'U',
    'U' = 'A',
    'G' = 'C',
    'C' = 'G',
    'Y' = 'R',
    'R' = 'Y',
    'S' = 'S',
    'W' = 'W',
    'K' = 'M',
    'M' = 'K',
    'B' = 'V',
    'D' = 'H',
    'H' = 'D',
    'V' = 'B',
    'N' = 'N',
    'a' = 'u',
    'u' = 'a',
    'g' = 'c',
    'c' = 'g',
    'y' = 'r',
    'r' = 'y',
    's' = 's',
    'w' = 'w',
    'k' = 'm',
    'm' = 'k',
    'b' = 'v',
    'd' = 'h',
    'h' = 'd',
    'v' = 'b',
    'n' = 'n'
  )
