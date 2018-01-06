#' Translate DNA or RNA
#'
#' Translate a vector of nucleotide sequences into a vector of protein
#' sequences.
#'
#' @param x `[character]`
#'
#'   Sequences to translate. Must only contain valid IUPAC DNA or RNA
#'   characters.
#'
#' @param type `[function]`
#'
#'   Used to coerce x to a class compatible with [Biostrings::translate].
#'   Typically one of [Biostrings::DNAStringSet], or
#'   [Biostrings::RNAStringSet]. Defaults to [Biostrings::DNAStringSet].
#'
#' @param translation `[named:character(64)]`
#'
#'   Amino acid letters with names corresponding to each of the 64 possible DNA
#'   codons. Defaults to `TRANSLATION`.
#'
#' @param start `[character]`
#'
#'   Start codons. The first occurence of one of these codons will be
#'   translated to `'M'`, and further occurences will be translated to that
#'   defined by `translation`. Defaults to `'ATG'`.
#'
#' @param abiguities `['solve'|'error'|'error.if.X'|'X']`
#'
#'   How ambiguities be handled. See [Biostrings::translate] for more details.
#'   Defaults to `'solve'` which will translate to the appropriate amino acid,
#'   and will translate to `'X'` if not possible to unambiguously translate to
#'   a single amino acid.
#'
#' @details
#'   This function is a simple wrapper around [Biostrings::translate] that
#'   takes a character vector as input and a returns a character vector.
#'   The default behavior is slightly different in that ambiguous bases will
#'   be solved rather than causing an error, and alternative start codons are
#'   not used without being explicitly specified with the `start` parameter.
#'   This prevents a common "gotcha" when translating a subsequence with
#'   [Biostrings::translate] using the default [Biostrings::GENETIC_CODE]
#'   where it would unexpectedly translate `"TTG"`, or `"CTG"` to `"M"` rather
#'   than `"L"` without any warning. To reproduce this behavior, simply specify
#'   `start = c('ATG', 'TTG', 'CTG')` which is only sensible if you are
#'   translating an entire coding sequence.
#'
#'
#' @return Returns a character vector with the same length as `x`.
#'
#' @examples
#'   # Default settings assume DNA input and will attempt to solve ambiguity
#'   translate(c('ATGAAAGGRTGA', 'ATGGGRRAATAA'))
#'
#'   # RNA can also be translated
#'   translate('AUGAAAGGRUGA', type = Biostrings::RNAStringSet)
#'
#'   # Alternative start codons may be specified
#'   translate(c('AAAATGAAA', 'ATGAAAAAA'), start = c('ATG', 'AAA'))
#'
#'   # Alternative translation code can be specified
#'   codon_code <- TRANSLATION
#'   codon_code['TGG'] <- 'R'
#'   translate(c('AAATGGTGG', 'ATGTGGTGG'), translation = codon_code)
#'
#'   # Other common translation tables
#'   Biostrings::GENETIC_CODE_TABLE
#'   Biostrings::getGeneticCode('Yeast Mitochondrial', full.search = TRUE)
#'
#' @importFrom Biostrings DNAStringSet RNAStringSet
#' @export
#' @md

translate <- function(x,
                      type = Biostrings::DNAStringSet,
                      translation = TRANSLATION,
                      start = "ATG",
                      ambiguities = "solve") {

  attr(translation, "alt_init_codons") <- start
  aa <- Biostrings::translate(type(x), translation, ambiguities)
  return(as.character(aa))
}

#' @format **TRANSLATION** - A named `character` vector with length 64. Names
#' correspond to codons, values correspond to single letter amino acid codes,
#' and '*' for stop.
#'
#' @rdname translate
#' @export
#' @md

TRANSLATION <-
  c(
    TTT = "F", TTC = "F", TTA = "L", TTG = "L", TCT = "S", TCC = "S",
    TCA = "S", TCG = "S", TAT = "Y", TAC = "Y", TAA = "*", TAG = "*",
    TGT = "C", TGC = "C", TGA = "*", TGG = "W", CTT = "L", CTC = "L",
    CTA = "L", CTG = "L", CCT = "P", CCC = "P", CCA = "P", CCG = "P",
    CAT = "H", CAC = "H", CAA = "Q", CAG = "Q", CGT = "R", CGC = "R",
    CGA = "R", CGG = "R", ATT = "I", ATC = "I", ATA = "I", ATG = "M",
    ACT = "T", ACC = "T", ACA = "T", ACG = "T", AAT = "N", AAC = "N",
    AAA = "K", AAG = "K", AGT = "S", AGC = "S", AGA = "R", AGG = "R",
    GTT = "V", GTC = "V", GTA = "V", GTG = "V", GCT = "A", GCC = "A",
    GCA = "A", GCG = "A", GAT = "D", GAC = "D", GAA = "E", GAG = "E",
    GGT = "G", GGC = "G", GGA = "G", GGG = "G"
  )
