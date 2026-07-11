#' @keywords internal
.trapz <- function(x, y, check = TRUE) {
  if (check) {
    if (length(x) != length(y)) stop("x and y must have the same length.")
    if (length(x) < 2L) stop("Need at least two points for trapezoidal rule.")
    if (anyNA(x) || anyNA(y)) stop("x and y must not contain NA.")
    # optional: ensure increasing grid
    if (is.unsorted(x)) {
      o <- order(x)
      x <- x[o]; y <- y[o]
    }
  }
  sum((y[-1] + y[-length(y)]) * diff(x) / 2)
}
