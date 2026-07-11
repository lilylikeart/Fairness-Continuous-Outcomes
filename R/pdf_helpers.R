#' @importFrom stats quantile sd qnorm pnorm setNames
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @keywords internal
.pointwise_hellinger <- function(f0, f1) {
  # Hellinger integrand per point: |sqrt(f0) - sqrt(f1)| / sqrt(2)
  abs(sqrt(pmax(f0, 0)) - sqrt(pmax(f1, 0))) / sqrt(2)
}

#' @keywords internal
.kde_align_to_grid <- function(y0, y1, n = 1000,
                               align = "approx", adjustBW=1,
                               from = NULL, to = NULL) {
  align <- match.arg(align)

  if (align == "density") {
    # pick support if not supplied
    if (is.null(from) || is.null(to)) {
      from <- min(c(y0, y1), na.rm = TRUE)
      to   <- max(c(y0, y1), na.rm = TRUE)
    }
    if (!is.finite(from) || !is.finite(to) || from >= to) {
      stop("Invalid 'from'/'to' for KDE support.")
    }

    # both densities directly on identical grid
    d0 <- stats::density(y0, from = from, to = to, n = n)
    d1 <- stats::density(y1, from = from, to = to, n = n)

    xs <- d0$x
    f0 <- d0$y
    f1 <- d1$y

  } else {

  # KDEs on each group (defaults: Gaussian kernel, bw.nrd0)
  d0 <- stats::density(y0,adjust=adjustBW)
  d1 <- stats::density(y1,adjust=adjustBW)

  # Common grid spanning both supports
  xs <- seq(min(c(d0$x, d1$x)), max(c(d0$x, d1$x)), length.out = n)
  f0 <- stats::approx(d0$x, d0$y, xout = xs, rule = 2)$y
  f1 <- stats::approx(d1$x, d1$y, xout = xs, rule = 2)$y

  # for consistency, set from/to to the grid endpoints we actually used
  from <- xs[1]
  to   <- xs[length(xs)]
  }

  list(
    xs = xs, f0 = f0, f1 = f1,          # aligned grid and densities
    d0 = d0, d1 = d1,                   # original density() objects
    align = align,                      # which method was used
    from = from, to = to                # support endpoints used
  )
}

#' @keywords internal
.density_fun_from_kde <- function(d) {
  # Returns a function f(x) based on linear interpolation of a KDE
  function(x) {
    stats::approx(d$x, d$y, xout = x, rule = 2)$y
  }
}

#' @keywords internal
.slice_window <- function(xs, a, b) {
  which(xs >= a & xs <= b)
}


# internal
.bin_numeric <- function(x, method = c("hist","boxplot"), num_bins = 10) {
  method <- match.arg(method)
  x <- x[!is.na(x)]
  if (length(x) < 2L) stop("Not enough non-NA values to bin.")
  if (method == "hist") {
    r <- range(x); breaks <- seq(r[1], r[2], length.out = num_bins + 1L)
  } else {
    #qs <- stats::quantile(x, c(.25,.5,.75), na.rm = TRUE)
    #qs <- boxplot_stats(x)
    #breaks <- c(-Inf, qs, Inf)
    qs <- stats::quantile(
      x,
      probs = seq(0, 1, length.out = num_bins + 1),
      na.rm = TRUE
    )
    # remove the 0 and 1 boundaries (since you'll add -Inf / Inf manually)
    qs <- qs[-c(1, length(qs))]

    breaks <- c(-Inf, qs, Inf)
  }
  list(
    breaks = breaks,
    bins = cut(x, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  )
}



