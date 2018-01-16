#' @importFrom magrittr %>%
#' @importFrom dplyr group_by arrange mutate ungroup
#' @export

add_exon_details <- function(cds) {
  next_frame <- c(2, 3, 1)

  cds %>%
    group_by(tx) %>%
    arrange(tx, exon) %>%
    mutate(
      exon_width    = (end - start) + 1L,
      exon_position = cumsum(lag(exon_width, default = 0L)) + 1L,
      exon_frame    = as.integer(exon_position %% 3) + 1L
    ) %>%
    ungroup()
}

#' @importFrom tidyr separate_rows
#' @importFrom fuzzyjoin genome_inner_join
#' @importFrom dplyr select
#' @export

inner_join_cds_vcf <- function(cds,
                               vcf,
                               cds_keys = c('chromosome' = 'chr'   , 'start' = 'start', 'end' = 'end'),
                               vcf_keys = c('chromosome' = 'CHROM' , 'start' = 'POS',   'end' = 'POS'),
                               separate_alleles = TRUE) {

  if (separate_alleles) vcf <- tidyr::separate_rows(vcf, 'ALT', sep = ',')

  cds$chromosome <- cds[[cds_keys['chromosome']]]
  cds$start      <- cds[[cds_keys['start']]]
  cds$end        <- cds[[cds_keys['end']]]
  vcf$chromosome <- vcf[[vcf_keys['chromosome']]]
  vcf$start      <- vcf[[vcf_keys['start']]]
  vcf$end        <- vcf[[vcf_keys['end']]]

  if (!any(unique(cds$chromosome) %in% unique(vcf$chromosome))) {
    stop('There were no matching chromosome names between cds and vcf tables.')
  }

  fuzzyjoin::genome_inner_join(cds, vcf, by = c('chromosome', 'start', 'end')) %>%
    select(
      # CDS columns
      tx, exon, chr, exon_strand = strand,
      exon_start = start.x, exon_end = end.x,
      vcf_start = start.y, vcf_end = end.y,
      # Typical VCF columns
      ID, CHROM, POS, REF, ALT, QUAL, FILTER, INFO,
      everything(),
      -chromosome.x, -chromosome.y
    )
}


test <- function() {
  genome <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens

  vcf <-
    read_vcf('~/GitServer/Manuscript-ABE/data/ClinVar/clinvar.vcf.gz') %>%
    mutate(CHROM = paste0('chr', gsub('MT', 'M', CHROM)))

  cds <-
    read_csv('~/GitServer/Manuscript-ABE/data/CDS/Hsapiens-UCSC-hg38.csv') %>%
    filter(chr %in% paste0('chr', c(1:22, 'X', 'Y', 'M'))) %>%
    select(chr, start, end, strand, tx, gene, exon) %>%
    add_exon_details()

  vcf_with_exons <- inner_join_cds_vcf(cds, vcf)
}
