#' Add variant effect prediction
#'
#' Experimental, use with caution
#'
#' @param cds Coding sequence coordinates
#' @param vcf Data frame of variants from a variant call format file
#' @param genome A genome reference compatible with [get_genomic_sequence]
#'
#' @details
#'
#' `cds`, `vcf`, and `genome` should have matching chromosome naming
#' conventions.
#'
#' `cds` must contain the following columns (rows represent exon coordinates):
#' `tx`, `exon`, `chr`, `strand`, `start`, and `end`
#'
#' `vcf` must contain the following standard VCF columns (rows represent variants):
#' `CHROM`, `POS`, `REF`, `ALT`.
#'
#' @examples
#' if (requireNamespace('BSgenome.Hsapiens.UCSC.hg38')) {
#'   library(tidyverse)
#'   library(mutagenesis)
#'
#'   # Example files provided in this package
#'   cds_file <- system.file('extdata/CDS.csv', package = 'mutagenesis')
#'   vcf_file <- system.file('extdata/VCF.vcf', package = 'mutagenesis')
#'
#'   # Coding sequence coordinates, variants, and a reference genome
#'   # are required to predict variant effects
#'   cds    <- read_csv(cds_file)
#'   vcf    <- read_vcf(vcf_file)
#'   genome <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
#'   vep    <- predict_variant_effect(cds, vcf, genome)
#'
#'   # An example summary of effects
#'   vep %>%
#'     select(
#'       gene, ID:INFO, ref_cds, alt_cds, ref_aa, alt_aa, mutation_type,
#'       exon_boundary_dist, vcf_start, vcf_end, exon_start, exon_end
#'     ) %>%
#'     distinct() %>%
#'     count(ref_aa, mutation_type) %>%
#'     arrange(mutation_type, -n)
#' }
#'
#' @export
#' @md

predict_variant_effect <- function(cds, vcf, genome) {

  # Prepare cds and vcf tables
  if (!hasName(cds, 'exon_frame')) cds <- add_exon_details(cds)

  vcf <-
    vcf %>%
    mutate(
      # Used to track variants after join
      vcf_id      = 1:n(),
      # Width of ref used to determine CDS join (i.e. start/end)
      ref_length  = str_length(REF),
      alt_length  = str_length(ALT),
      # Define VCF end and start based on ref_length
      # Caveats: fails for complex variants and missing should be encoded with
      # an empty string.
      start       = POS,
      end         = start + (ref_length - 1L)
    )

  # Join vcf and cds exons by overlapping coordinates
  cds_keys <- c('chromosome' = 'chr',    'start' = 'start', 'end' = 'end')
  vcf_keys <- c('chromosome' = 'CHROM' , 'start' = 'start', 'end' = 'end')
  joined <- mutagenesis::inner_join_cds_vcf(cds, vcf, cds_keys, vcf_keys)

  # Annotate variants
  joined %>%
    mutate(
      # Determine distance of Ref and alt boundaries to 5' of exon
      exon_5prime       = ifelse(exon_strand == '+', exon_start, exon_end),
      vcf_5prime        = ifelse(exon_strand == '+', vcf_start, vcf_end),
      vcf_3prime        = ifelse(exon_strand == '+', vcf_end, vcf_start),
      alt_5prime        = vcf_5prime,
      # used only to calculate frame of 3' end relative to exon 5' end
      alt_3prime        = ifelse(exon_strand == '+', vcf_start + (alt_length - 1L), vcf_end - (alt_length - 1L)),
      # Given distance to exon's 5' end determine frame for REF and ALT boundaries
      ref_5prime_frame  = get_frame(abs(vcf_5prime - exon_5prime) + (exon_frame - 1L)),
      ref_3prime_frame  = get_frame(abs(vcf_3prime - exon_5prime) + (exon_frame - 1L)),
      alt_5prime_frame  = ref_5prime_frame,
      alt_3prime_frame  = get_frame(abs(alt_3prime - exon_5prime) + (exon_frame - 1L)),
      # Given frame of boundaries, determine range for REF and ALT with complete codons
      inframe_ref_start = ifelse(exon_strand == '+', vcf_start + c(0, -1, -2)[ref_5prime_frame], vcf_start + c(-2, -1, 0)[ref_3prime_frame]),
      inframe_ref_end   = ifelse(exon_strand == '+', vcf_end   + c(2,  1,  0)[ref_3prime_frame], vcf_end   + c( 0,  1, 2)[ref_5prime_frame]),
      inframe_alt_start = ifelse(exon_strand == '+', vcf_start + c(0, -1, -2)[alt_5prime_frame], vcf_start + c(-2, -1, 0)[alt_3prime_frame]),
      inframe_alt_end   = ifelse(exon_strand == '+', vcf_end   + c(2,  1,  0)[alt_3prime_frame], vcf_end   + c( 0,  1, 2)[alt_5prime_frame]),
      # Extract boundary genomic sequence for REF and ALT
      ref_up = ifelse(vcf_start - 1L < inframe_ref_start, '', get_genomic_sequence(CHROM, '+', inframe_ref_start, vcf_start - 1L,  genome)),
      ref_dn = ifelse(vcf_end   + 1L > inframe_ref_end,   '', get_genomic_sequence(CHROM, '+', vcf_end   + 1L,    inframe_ref_end, genome)),
      alt_up = ifelse(vcf_start - 1L < inframe_alt_start, '', get_genomic_sequence(CHROM, '+', inframe_alt_start, vcf_start - 1L,  genome)),
      alt_dn = ifelse(vcf_end   + 1L > inframe_alt_end,   '', get_genomic_sequence(CHROM, '+', vcf_end   + 1L,    inframe_alt_end, genome)),
      ref_seq = str_c(tolower(ref_up), REF, tolower(ref_dn)),
      alt_seq = str_c(tolower(alt_up), ALT, tolower(alt_dn)),

      # Variant types describe the general change in genomic sequence
      variant_type = case_when(
        str_detect(ALT, '[<>\\[\\]]')       ~ 'Complex',
        ref_length == 1L & alt_length == 1L ~ 'Single nucleotide',
        ref_length == 2L & alt_length == 2L ~ 'Dinucleotide',
        ref_length == 3L & alt_length == 3L ~ 'Trinucleotide',
        ref_length == alt_length            ~ 'Substitution',
        ref_length != alt_length            ~ 'Indel'
      ),
      # Minimum distance to exon boundary, -1 if crosses junction, missing if complex
      exon_boundary_dist = case_when(
        # No attempt to determine distance for complex variants
        variant_type == 'Complex' ~ NA_integer_,
        vcf_start < exon_start | vcf_end > exon_end ~ -1L, # overlaps a junction
        TRUE ~ pmin(
          abs(vcf_start - exon_start),
          abs(vcf_start - exon_end),
          abs(vcf_end   - exon_start),
          abs(vcf_end   - exon_end)
        )
      ),
      splicing_type = case_when(
        variant_type == 'Complex' ~ NA_character_,
        exon_boundary_dist < 0L   ~ 'Overlaps splice junction',
        exon_boundary_dist < 2L   ~ 'Adjacent to splice junction',
        TRUE                      ~ 'Within exon'
      ),
      ref_cds = ifelse(exon_strand == '+', ref_seq, reverse_complement(ref_seq)),
      alt_cds = ifelse(exon_strand == '+', alt_seq, reverse_complement(alt_seq)),
      ref_aa  = ifelse(splicing_type == 'Within exon', mutagenesis::translate(ifelse(splicing_type == 'Within exon', ref_cds, '')), NA_character_),
      alt_aa  = ifelse(splicing_type == 'Within exon', mutagenesis::translate(ifelse(splicing_type == 'Within exon', alt_cds, '')), NA_character_),
      mutation_type = case_when(
        variant_type == 'Complex'                                     ~ 'Complex',
        splicing_type != 'Within exon'                                ~ 'Splicing',
        variant_type == 'Indel' & (ref_length - alt_length) %% 3 != 0 ~ 'Frameshift',
        ref_aa == alt_aa                                              ~ 'Silent',
        ref_aa != alt_aa & str_detect(alt_aa, '[*]')                  ~ 'Nonsense',
        ref_aa != alt_aa                                              ~ 'Missense'
      )
    ) %>%
    select(
      -exon_5prime, -vcf_5prime, -vcf_3prime, -alt_5prime, -alt_3prime,
      -ref_5prime_frame, -ref_3prime_frame, -alt_5prime_frame, -alt_3prime_frame,
      -ref_up, -ref_dn, -alt_up, -alt_dn
    )
}

get_frame <- function(x) {
  frame <- c(3, 1, 2)  # Multiple of 3 is frame 3 followed by frame 1 then 2
  frame[(x %% 3) + 1L] # 1 is added to allow indexing of `frame`
}
