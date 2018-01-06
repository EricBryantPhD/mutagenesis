#' Identify coding variants in a VCF file.
#'
#'
#'
#' @importFrom tidyr separate_rows
#' @importFrom magrittr %>%
#' @export
#' @md

coding_variants <- function(vcf, cds, genome, separate_alleles = TRUE) {

  # Extract header information to add back later
  vcf_header <- attr(vcf, 'vcf')

  # Sometimes multiple alleles are stored in the ALT field separated by commas
  if (separate_alleles) {
    vcf <- tidyr::separate_rows(vcf, 'ALT', sep = ',', convert = FALSE)
  }

  # Extract currently unsupported variant types and provide a warning
  unsupported <- filter(vcf, stringr::str_detect(ALT, '[^ATGCatgc]'))

  if (nrow(unsupported)) {
    warning(
      nrow(unsupported),
      'unsupported variants have been ignored (i.e. ambiguous sequences and structural variants). ',
      'Use attr(x, "vcf.unsupported") to access these variants'
    )
  }

  vcf <-
    setdiff(vcf, unsupported) %>%
    mutate(start = POS, end = POS + (stringr::str_length(REF) - 1L))

  cds <- add_exon_frame(cds)

  # Annotate positions with affected CDS transcripts
  coding <-
    fuzzyjoin::genome_inner_join(vcf, cds, by = c('CHROM', 'start', 'end')) %>%
    select(
      ID, CHROM = CHROM.x,
      REF_start = start.x, REF_end = end.x, REF, ALT,
      QUAL, FILTER, INFO,
      tx, exon,
      exon_start = start.y, exon_end = end.y, exon_strand = strand, exon_frame
    ) %>%
    mutate(REF_check = get_sequence_range(genome, CHROM, strand = '+', REF_start, REF_end))

  test <-
    coding %>%
    group_by(tx) %>%
    arrange(tx, exon) %>%
    summarise(
      REF_allele = get_sequence_range(genome, CHROM, exon_strand, REF_start, REF_end),
      ALT_allele = get_sequence_range(genome, CHROM, exon_strand, REF_start, REF_end)
    )

  # Reset and create attributes
  attr(coding, 'vcf.unsupported') <- unsupported
  attr(coding, 'vcf') <- vcf_header
  class(coding) <- c('vcf', class(coding))
  return(coding)
}



neighboring_coding_sequence <- function(txs, REF_start, REF_end, cds) {
  cds %>%
    filter(tx %in% txs, at >= start, at <= end) %>%
    group_by(tx) %>%
    arrange(tx, exon) %>%
    summarise(
      starts = str_c(start, ','),
      ends = str_c(end, ','),
      coding_sequence = get_sequence_range(genome, CHROM, strand, start, end)
    )
}

add_exon_frame <- function(cds) {
  next_frame <- c(2, 3, 1)

  cds %>%
    group_by(tx) %>%
    arrange(tx, exon) %>%
    mutate(
      width      = (end - start) + 1L,
      n_upto     = cumsum(lag(width, default = 0L)),
      exon_frame = as.integer(n_upto %% 3),
      exon_frame = next_frame[ifelse(exon_frame == 0L, 3L, exon_frame)]
    ) %>%
    ungroup() %>%
    select(-width, -n_upto)
}

coding_frame <- function(exon_frame)

test <- function() {
  genome <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens

  vcf <-
    read_vcf('~/GitServer/Manuscript-ABE/data/ClinVar/clinvar.vcf.gz') %>%
    mutate(CHROM = paste0('chr', gsub('MT', 'M', CHROM)))

  cds <-
    read_csv('~/GitServer/Manuscript-ABE/data/CDS/Hsapiens-UCSC-hg38.csv') %>%
    filter(chr %in% paste0('chr', c(1:22, 'X', 'Y', 'M'))) %>%
    select(CHROM = chr, start, end, strand, tx, gene, exon)

  coding_variants(vcf, cds, genome)
}
