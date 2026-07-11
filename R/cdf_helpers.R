#' @noRd
.weights_or_one <- function(n, w) {
  if (is.null(w)) return(rep(1, n))
  if (length(w) != n) w <- rep(w, length.out = n)
  w[is.na(w) | !is.finite(w)] <- 0
  w
}

# Weighted ECDF as a right-continuous step function (approxfun) adapted from stat-ecdf.R in ggplot
#' @noRd
.ecdf<- function(x, weights = NULL) {

  if (length(x) == 0L) stop("Empty vector for ECDF.")

  weights <- weights %||% 1
  weights <- vctrs::vec_recycle(weights, length(x))

  if (!all(is.finite(weights))) {
    cli::cli_warn(c(paste0(
      "The {.field weight} aesthetic does not support non-finite or ",
      "{.code NA} values."
    ), "i" = "These weights were replaced by {.val 0}."))
    weights[!is.finite(weights)] <- 0
  }

  ok <- !(is.na(x) | is.na(weights))
  x <- x[ok]
  weights <- weights[ok]

  ord <- order(x)
  x <- x[ord]
  weights <- weights[ord]

  total <- sum(weights)

  if (abs(total) < 1000 * .Machine$double.eps) {
    if (total == 0) {
      cli::cli_abort(paste0(
        "Cannot compute eCDF when the {.field weight} aesthetic sums up to ",
        "{.val 0}."
      ))
    }
    cli::cli_warn(c(
      "The sum of the {.field weight} aesthetic is close to {.val 0}.",
      "i" = "Computed eCDF might be unstable."
    ))
  }

  ux <- unique(x)
  idx <- match(x, ux)
  #wsum <- tapply(w, idx, sum)
  #wsum <- rowsum(weights, idx, reorder = FALSE)
  #wsum <- as.numeric(wsum)  # keep order 1..k

  wsum <- vapply(
    split(weights, idx),
    sum, numeric(1)
  )

  csum <- cumsum(wsum) / total

  stats::approxfun(ux, csum,
                   method = "constant",
                   yleft = 0, yright = 1,
                   f = 0, ties = "ordered")
}

# Build evaluation grid for CDFs
# If n is NULL: exact knot grid (union of unique values) → exact integral for step CDFs
# Else: uniform grid over observed range
#' @noRd
.make_cdf_grid <- function(y0, y1, n = NULL) {
  if (is.null(n)) {
    sort(unique(c(y0, y1)))
  } else {
    rng <- range(c(y0, y1), na.rm = TRUE)
    seq(rng[1], rng[2], length.out = n)
  }
}
