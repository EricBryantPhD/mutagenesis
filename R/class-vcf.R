#' @export
print.vcf <- function(x) {
  cat(attr(x, 'vcf'), sep = '\n')
  cat('\n')
  tibble:::print.tbl(x)
}

#' @export
filter.vcf <- function(.data, ...) reclass(.data, NextMethod())

# Reclass maintains class and attributes of an object after calling
# `NextMethod()`
reclass <- function(x, result) UseMethod('reclass')

reclass.default <- function(x, result) {
  class(result) <- unique(c(class(x)[1], class(result)))
  attr(result, class(x)[1]) <- attr(x, class(x)[1])
  result
}
