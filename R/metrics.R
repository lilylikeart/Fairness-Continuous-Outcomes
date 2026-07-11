#' Areas between PDFs + Hellinger (with optional local window around 0)
#'
#' Computes:
#' - AreaPDF (L1 distance between group PDFs) via trapezoidal integration
#' - Hellinger distance (global)
#' - Optional probability mass for each group in a tiny window around 0
#' - Pointwise Hellinger values (useful for plotting)
#'
#' @param data A data.frame.
#' @param target_col Character; numeric target column name.
#' @param group_col Character; binary group column name (exactly two groups).
#' @param n Integer; number of grid points for the common support (default 1000).
#' @param align character, "density" or "approx"
#' @param bwAdj Bandwidth adjustment factor passed to
#'   \code{\link[stats]{density}}.
#' @param from,to numeric, optional explicit support when align = "density"
#' @param verbose Logical; if TRUE, prints brief sizes/metrics (default FALSE).
#' @return A list with metrics and aligned densities.
#' @export
compute_AreaPDF <- function(data, target_col, group_col, n = 1000,
                             align = "approx", bwAdj=1,
                             from = NULL, to = NULL,
                             verbose = FALSE) {

  stopifnot(target_col %in% names(data), group_col %in% names(data))

  align <- match.arg(align, c("density", "approx"))

  target <- data[[target_col]]
  group  <- data[[group_col]]
  gs <- sort(unique(group))
  if (!is.numeric(group) || length(gs) != 2L || !setequal(gs, c(0, 1))) {
    stop("group_col must contain exactly two levels coded as numeric 0 and 1.")
  }
  y0 <- target[group == 0]
  y1 <- target[group == 1]

  if (verbose) {
    message(sprintf("n_%s = %d, n_%s = %d",
                    as.character(gs[1]), length(y0),
                    as.character(gs[2]), length(y1)))
  }

  # KDEs aligned onto a common grid
  kd <- .kde_align_to_grid(y0, y1, n = n, align, bwAdj, from, to)
  xs <- kd$xs; f0 <- kd$f0; f1 <- kd$f1

  # L1 area between PDFs (AreaPDF)
  area_pdf <- .trapz(xs, abs(f0 - f1))

  #area_pdf<-sum(abs(f0 -f1)* mean(diff(xs)))

  # Global Hellinger distance (integral of pointwise term)
  h_vals <- .pointwise_hellinger(f0, f1)
  hellinger <- sqrt(.trapz(xs, (sqrt(pmax(f0,0)) - sqrt(pmax(f1,0)))^2)) / sqrt(2)

  # Optional: probability mass in a tiny window around 0 for each group
  eps_around0 <- 1e-3
  f0_fun <- .density_fun_from_kde(kd$d0)
  f1_fun <- .density_fun_from_kde(kd$d1)
  prob0_around0 <- stats::integrate(f0_fun, lower = -eps_around0, upper = eps_around0)$value
  prob1_around0 <- stats::integrate(f1_fun, lower = -eps_around0, upper = eps_around0)$value

  # Also compute Hellinger restricted to the small window (useful diagnostic)
  idx_win <- .slice_window(xs, -eps_around0, eps_around0)
  h_local <- if (length(idx_win) >= 2) {
    .trapz(xs[idx_win], .pointwise_hellinger(f0[idx_win], f1[idx_win]))
  } else NA_real_

  if (verbose) {
    message(sprintf("AreaPDF = %.6f; Hellinger = %.6f;
                    prob0_around0=%.6f, prob1_around0=%.6f; Hellinger_around0 =%.6f",
                    area_pdf, hellinger, prob0_around0, prob1_around0, h_local))
  }

  list(
    AreaPDF = area_pdf,
    Hellinger = hellinger,
    Hellinger_window = h_local,
    grid = xs,
    density0 = f0,
    density1 = f1,
    pointwise_hellinger = h_vals,
    prob_at_0_group0 = prob0_around0,
    prob_at_0_group1 = prob1_around0
  )
}

#' Area between CDFs (AreaCDF) for two groups
#'
#' Computes the L1 area between two empirical CDFs. In 1D this equals the
#' first Wasserstein distance W1. Supports optional per-row weights.
#'
#' @param data A data.frame.
#' @param target_col Character; numeric target column.
#' @param group_col Character; binary group column (two groups; values can be any two labels).
#' @param n Integer or NULL; grid size. If NULL, uses exact knots (union of unique values).
#' @param weights_col Optional character; column name of nonnegative weights.
#' @param verbose Logical; print small diagnostics.
#' @return A list with \code{AreaCDF}, the evaluation \code{grid}, and the two CDF
#'   vectors on that grid (named by the original group labels).
#' @examples
#' set.seed(1)
#' df <- data.frame(score = c(rnorm(150,0), rnorm(150,0.5)),
#'                  grp   = rep(c(0,1), each = 150))
#' out <- compute_AreaCDF(df, "score", "grp", n = NULL)
#' out$AreaCDF
#' @export
compute_AreaCDF <- function(data, target_col, group_col,
                             n = NULL, weights_col = NULL,
                             verbose = FALSE) {
  stopifnot(target_col %in% names(data), group_col %in% names(data))

  target <- data[[target_col]]
  group <- data[[group_col]]

  # drop NAs
  ok <- !(is.na(target) | is.na(group))
  if (!all(ok)) {
    if (verbose) message(sprintf("Dropping %d rows with NA.", sum(!ok)))
    target <- target[ok]; group <- group[ok]
    if (!is.null(weights_col)) data <- data[ok, , drop = FALSE]
  }

  # resolve group labels (allow any two labels)
  #glabs <- unique(as.character(group))
  #if (length(glabs) != 2L) stop("group_col must have exactly two groups.")
  if (!all(sort(unique(group)) == c(0,1)))
    stop("group_col must contain 0 and 1.")

  y0 <- target[group == 0]
  y1 <- target[group == 1]

  # optional weights per row
  w0 <- w1 <- NULL
  if (!is.null(weights_col)) {
    w <- data[[weights_col]]
    w0 <- w[group == 0]
    w1 <- w[group == 1]
  }

  # weighted ECDFs
  F0 <- .ecdf(y0, w0)
  F1 <- .ecdf(y1, w1)

  # evaluation grid
  xs <- .make_cdf_grid(y0, y1, n = n)

  c0 <- F0(xs)
  c1 <- F1(xs)

  # integral of |F0 - F1| over R; outside [min,max] both CDFs are equal (0 or 1)
  area_cdf <- .trapz(xs, abs(c0 - c1))

  #dx <- diff(xs)

  #print(paste("ecdf_values_g0:",ecdf_values_g0))
  #print(paste("ecdf_values_g1:",ecdf_values_g1))

  #absolute_differences <- abs(c0 - c1)

  #area_cdf <- (absolute_differences[-length(absolute_differences)] + absolute_differences[-1]) / 2 * dx

  #sum_area_cdf <- sum(area_cdf)

  if (verbose) {
    message(sprintf("area_cdf (W1) = %.6f, grid points = %d", area_cdf,length(xs)))
  }

  list(
    AreaCDF = area_cdf,
    grid = xs,
    ecdf0 = stats::setNames(c0, NULL),
    ecdf1 = stats::setNames(c1, NULL)
  )
}

#' Compute Conditional AreaPDF Fairness Metric
#'
#' Computes the AreaPDF fairness metric within bins of a conditioning
#' variable and summarizes the fairness across all bins.
#'
#' The conditioning variable is partitioned into bins using either
#' equal-width histogram bins or boxplot-style breaks. Within each bin,
#' the AreaPDF and Hellinger distance are computed between the two groups.
#' Bins with fewer than `min_per_group` observations in either group are
#' skipped.
#'
#' @param data A data frame containing the observations.
#' @param target_col Name of the numeric outcome (score) column.
#' @param condition_col Name of the conditioning variable used for binning.
#' @param group_col Name of the binary group column (0/1).
#' @param bwAdj Bandwidth adjustment passed to `compute_AreaPDF()`.
#' @param method Method used to construct bins. Either `"hist"` or
#'   `"boxplot"`.
#' @param num_bins Number of bins when `method = "hist"`.
#' @param min_per_group Minimum number of observations required in each
#'   group within a bin to compute the metric.
#' @param n Number of grid points used for density estimation.
#' @param verbose Logical; if `TRUE`, progress information is printed.
#'
#' @return A list containing:
#' \describe{
#'   \item{mean_AreaPDF}{Mean AreaPDF across valid bins.}
#'   \item{breaks}{Bin boundaries.}
#'   \item{per_bin}{Summary statistics for each bin.}
#'   \item{pointwise}{Pointwise Hellinger distance and densities.}
#'   \item{used_data}{Input data annotated with bin assignments.}
#'   \item{groups}{Group labels.}
#' }
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   score = rnorm(200),
#'   age = runif(200, 18, 65),
#'   group = rep(c(0, 1), each = 100)
#' )
#'
#' compute_AreaPDF_cond(
#'   df,
#'   target_col = "score",
#'   condition_col = "age",
#'   group_col = "group"
#' )
#'
#' @export
compute_AreaPDF_cond <- function(data, target_col, condition_col, group_col,bwAdj=1,
                                    method = c("hist","boxplot"),
                                    num_bins = 10,
                                    min_per_group = 2,
                                    n = 1000,
                                    verbose = FALSE) {
  method <- match.arg(method)

  # Extract vectors (keep NA rows; we'll mask per-bin)
  y      <- data[[target_col]]
  cond   <- data[[condition_col]]
  group  <- data[[group_col]]

  # Build bins (on non-NA cond)
  binfo <- .bin_numeric(cond, method = method, num_bins = num_bins)
  # Recompute bins on full vector, so length matches data
  bins <- cut(cond, breaks = binfo$breaks, include.lowest = TRUE, labels = FALSE)

  # Iterate bins
  bin_ids <- sort(unique(bins[!is.na(bins)]))
  per_bin <- vector("list", length(bin_ids))
  pointwise <- vector("list", length(bin_ids))

  k <- 0L
  for (b in bin_ids) {
    sel <- which(bins == b & !is.na(y) & !is.na(group))

    if (!length(sel)) next
    sub <- data[sel, , drop = FALSE]

    n_b <- length(sel)
    n0  <- sum(group[sel] == 0)
    n1  <- sum(group[sel] == 1)

    # default: no metrics (too few per group or will be filled below)
    area_pdf  <- NA_real_
    hellinger <- NA_real_
    pw_df     <- NULL

    can_compute <- (n0 >= min_per_group && n1 >= min_per_group)

    if (can_compute) {
      abp <- compute_AreaPDF(sub, target_col, group_col, align = "approx",bwAdj=bwAdj)
      area_pdf  <- abp$AreaPDF
      hellinger <- abp$Hellinger

      if (verbose) {
        message(sprintf(
          "bin %s: n=%d (n0=%d, n1=%d), AreaPDF=%.6f",
          b, n_b, n0, n1, area_pdf
        ))
      }

      pw_df <- data.frame(
        bin = b,
        x   = abp$grid,
        h   = abp$pointwise_hellinger,
        f0  = abp$density0,
        f1  = abp$density1,
        stringsAsFactors = FALSE
      )
    } else if (verbose) {
      message(sprintf(
        "bin %s: n=%d (n0=%d, n1=%d) - insufficient per-group samples (min_per_group=%d); metrics set to NA",
        b, n_b, n0, n1, min_per_group
      ))
    }

    k <- k + 1L
    per_bin[[k]] <- data.frame(
      bin        = b,
      n_samples  = n_b,
      n_group0   = n0,
      n_group1   = n1,
      AreaPDF    = area_pdf,
      Hellinger  = hellinger,
      stringsAsFactors = FALSE
    )

    pointwise[[k]] <- pw_df
  }


    # Need both groups with at least min_per_group
    #if (n0 < min_per_group || n1 < min_per_group) next


    #abp <- compute_AreaPDF(sub, target_col, group_col, align="approx")
    #area_pdf  <- abp$AreaPDF
    #hellinger <- abp$Hellinger
    # expected names: AreaPDF, Hellinger, grid, pointwise_hellinger, density0, density1

    #if (verbose) message(sprintf("bin %s: n=%d (n0=%d, n1=%d), areaPDF=%.6f", b, n_b, n0, n1,abp$AreaPDF))

    # k <- k + 1L
    # per_bin[[k]] <- data.frame(
    #   bin = b,
    #   n_samples = n_b,
    #   n_group0  = n0,
    #   n_group1  = n1,
    #   AreaPDF   = area_pdf,
    #   Hellinger = hellinger,
    #   stringsAsFactors = FALSE
    # )
    #
    # pointwise[[k]] <- data.frame(
    #   bin  = b,
    #   x    = abp$grid,
    #   h    = abp$pointwise_hellinger,
    #   f0   = abp$density0,
    #   f1   = abp$density1,
    #   stringsAsFactors = FALSE
    # )
    # per_bin_df   <- if (k) do.call(rbind, per_bin[seq_len(k)]) else
    # data.frame(bin=integer(), n_samples=integer(), n_group0=integer(),
    #            n_group1=integer(), AreaPDF=numeric(), Hellinger=numeric())
    # pointwise_df <- if (k) do.call(rbind, pointwise[seq_len(k)]) else
    # data.frame(bin=integer(), x=numeric(), h=numeric(), f0=numeric(), f1=numeric())


  # collapse per_bin
  per_bin <- per_bin[seq_len(k)]
  per_bin_df <- if (k) {
    do.call(rbind, per_bin)
  } else {
    data.frame(
      bin       = integer(),
      n_samples = integer(),
      n_group0  = integer(),
      n_group1  = integer(),
      AreaPDF   = numeric(),
      Hellinger = numeric()
    )
  }

  # collapse pointwise (skip NULL entries)
  valid_pw <- pointwise[!vapply(pointwise, is.null, logical(1))]
  pointwise_df <- if (length(valid_pw)) {
    do.call(rbind, valid_pw)
  } else {
    data.frame(bin = integer(), x = numeric(), h = numeric(), f0 = numeric(), f1 = numeric())
  }

  # Simple aggregate (equal-weight mean across bins)
  #mean_AreaPDF <- mean(per_bin_df$AreaPDF)
  mean_AreaPDF <- if (nrow(per_bin_df)) {
    mean(per_bin_df$AreaPDF, na.rm = TRUE)
  } else {
    NA_real_
  }

  # Return the 'used data' view without dplyr
  used_data <- setNames(
    data.frame(
      y      = y,
      cond   = cond,
      group  = group,
      bin    = bins
    ),
    c(target_col, condition_col, group_col, paste0("bin_", condition_col))
  )

  if (verbose) {
    message(sprintf("Mean AreaPDF across bins = %.6f", mean_AreaPDF))
    message(sprintf("Breaks: %s", paste(signif(binfo$breaks, 5), collapse = ", ")))
  }

  list(
    mean_AreaPDF = mean_AreaPDF,
    breaks         = binfo$breaks,
    per_bin        = per_bin_df,
    pointwise      = pointwise_df,
    used_data      = used_data,
    groups         = group
  )
}
