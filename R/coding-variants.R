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
      tx, exon, exon_strand = strand,
      exon_start = start.x, exon_end = end.x,
      vcf_start  = start.y, vcf_end  = end.y,
      # Typical VCF columns
      CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO,
      everything(),
      # Chromosome recorded by VCF information
      -chromosome.x, -chromosome.y
    )
}
