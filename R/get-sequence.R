#' Get genomic sequences given ranges

#' @param chr `[character]`
#'
#'   Chromosome names. Must match names returned by `names(genome)`.
#'
#' @param strand `[character]`
#'
#'   Sequence strands (+|-).
#'
#' @param start `[integer]`
#'
#'   Start coordinates of ranges.
#'
#' @param end `[integer]`
#'
#'   End coordinates of ranges
#'
#' @param genome `[BSgenome|DNAStringSet]`
#'
#'   A reference genome. See Details.
#'
#' @details The reference genome can be either a `BSgenome` object from a
#' BSgenome reference package (see [BSgenome::BSgenome]), or a `DNAStringSet`
#' object (see [Biostrings::DNAStringSet]). `BSgenome` objects offer faster
#' sequence aquisition, but are limited to existing BSgenome packages
#' (see [BSgenome::available.genomes]) whereas `DNAStringSet` objects can
#' be easily created from any standard FASTA file using
#' [Biostrings::readDNAStringSet].
#'
#' @importFrom Biostrings getSeq
#' @importFrom tibble tibble
#' @export
#' @md

get_genomic_sequence <- function(chr, strand, start, end, genome) {
  tibble::tibble(seqnames = chr, strand, start, end) %>%
    as('GRanges') %>%
    Biostrings::getSeq(genome, .) %>%
    as.character()
}

#' @rdname get_genomic_sequence
#' @export

get_genomic_variant <- function(chr, strand, start, end, vcf, genome) {
  #end_min <- pmin(end, ref_end)
  #end_max <- pmax(end, ref_end)
  #start_min <- pmin(start, ref_start)
  #start_max <- pmax(start, ref_start)
  #alt_width <- stringr::str_length(alt) # TODO
  upstream <- get_sequence_range(genome, chr, '+', start, ref_start - 1L)
  dnstream <- get_sequence_range(genome, chr, '+', ref_end + 1L, end)
  plus_strand <- stringr::str_c(upstream, alt, dnstream)
  ifelse(strand == '+', plus_strand, reverse_complement(plus_strand))
}

#' @rdname get_genomic_sequence
#' @export

get_coding_sequence <- function(chr, strand, start, end, cds, genome) {

}

#' @rdname get_genomic_sequence
#' @export

get_coding_variant <- function(chr, strand, start, end, cds, vcf, genome) {

}
