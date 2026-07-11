#' Bootstrap Confidence Interval for Fairness Metrics
#'
#' Computes bootstrap confidence intervals for fairness metrics such as
#' AreaCDF, AreaPDF, and conditional AreaPDF.
#'
#' @param df A data frame containing the observations.
#' @param target Name of the target (score) column.
#' @param group Name of the binary group column.
#' @param cond Optional conditioning variable used when
#'   `metric = "AreaPDF_cond"`.
#' @param metric Fairness metric to compute. One of
#'   `"AreaCDF"`, `"AreaPDF"`, or `"AreaPDF_cond"`.
#' @param B Number of bootstrap resamples.
#' @param conf Confidence level.
#' @param seed Random seed.
#' @param perc_ci Logical; compute percentile confidence interval.
#' @param compute_bca Logical; compute BCa confidence interval.
#' @param studentized Logical; compute bootstrap-t confidence interval.
#' @param jackknife_max Maximum number of jackknife samples.
#'
#' @return A list containing the point estimate, bootstrap samples,
#' confidence intervals, and standard errors.
#'
#' @export
bootstrap_metric_ci <- function(df, target, group, cond = NULL, metric = c("AreaCDF","AreaPDF","AreaPDF_cond"),
                           B = 2000, conf = 0.95, seed = 1, perc_ci = TRUE, compute_bca = FALSE, studentized = TRUE,
                           jackknife_max = Inf) {
  set.seed(seed);
  metric <- match.arg(metric, choices = c("AreaCDF","AreaPDF","AreaPDF_cond"))

  # 1) Select the metric function once
  metric_fn <- switch(
    metric,
    AreaCDF      = function(data) compute_AreaCDF(data, target, group)$AreaCDF,
    AreaPDF      = function(data) compute_AreaPDF(data, target, group, align = "approx")$AreaPDF,
    AreaPDF_cond = function(data) compute_AreaPDF_cond(data, target, cond, group)$mean_AreaPDF
  )
  # 2) Point estimate
  est <- metric_fn(df)
  cat(sprintf("%s est: %.4f\n", metric, est))


  # 3) Stratified indices (assumes binary 0/1 group labels)
  g <- df[[group]];
  idx0 <- which(g == 0)
  idx1 <- which(g == 1)

  #est <- if (which == "AreaCDF") {compute_AreaCDF(df, target, group)$AreaCDF}
  #      else if (which == "AreaPDF"){compute_AreaPDF(df, target, group, align = "approx")$AreaPDF}
  #      else if (which == "AreaPDF_cond"){compute_AreaPDF_cond(df, target, group, align = "approx")$AreaPDF}
  #      else{stop("Invalid 'which' argument: must be one of 'AreaCDF', 'AreaPDF', or 'AreaKS'.")}
  out <- list(point = est)

  bootstrap <- numeric(B)
  if (perc_ci){
    for (b in seq_len(B)) {
      i <- c(sample(idx0, length(idx0), TRUE),
             sample(idx1, length(idx1), TRUE))
      d <- df[i, , drop = FALSE]
      bootstrap[b] <- metric_fn(d)

      #bootstrap[b] <- if (which == "AreaCDF") {compute_AreaCDF(df, target, group)$AreaCDF}
      #else if (which == "AreaPDF"){compute_AreaPDF(df, target, group, align = "approx")$AreaPDF}
      #else if (which == "AreaPDF_cond"){compute_AreaPDF_cond(df, target, group, align = "approx")$AreaPDF}
      #else{stop("Invalid 'which' argument: must be one of 'AreaCDF', 'AreaPDF', or 'AreaKS'.")}
    }
    alpha <- (1 - conf)/2
    ci_perc <- as.numeric(quantile(bootstrap, c(alpha, 1 - alpha), na.rm = TRUE))
    se_boot <- stats::sd(bootstrap, na.rm = TRUE)
    out$ci_perc <- ci_perc
    out$se <- se_boot
    out$bootstrap <- bootstrap
  }

  # --- optional: BCa CI ---
  if (compute_bca) {
    check_bootstrap_skew(bootstrap)
    # Bias-correction z0
    # proportion of draws < est (with 0.5 weight on ties), clamped away from 0/1
    R <- length(bootstrap)
    p0 <- (sum(bootstrap < est, na.rm = TRUE) + 0.5 * sum(bootstrap == est, na.rm = TRUE)) / R
    eps <- 1 / (R + 1)
    p0 <- min(max(p0, eps), 1 - eps)
    z0 <- stats::qnorm(p0)

    # Acceleration 'a' via leave-one-out jackknife
    n <- nrow(df)
    # Optionally subsample the jackknife for speed if n is huge
    jk_idx <- if (is.finite(jackknife_max) && jackknife_max < n) {
      sort(sample.int(n, jackknife_max))
    } else seq_len(n)

    t_j <- numeric(length(jk_idx))
    keep <- rep(TRUE, length(jk_idx))

    for (k in seq_along(jk_idx)) {
      ii <- jk_idx[k]
      d_loo <- df[-ii, , drop = FALSE]
      # guard: require both groups present; KDE needs >= 2 per group ideally
      gg <- d_loo[[group]]
      if (length(unique(gg)) < 2L || min(table(gg)) < 1L) { keep[k] <- FALSE; next }
      if (metric == "AreaCDF") {
        t_j[k] <- compute_AreaCDF(d_loo, target, group)$AreaCDF
      } else {
        # for KDE, insist on at least 2 per group; otherwise skip this jackknife point
        if (min(table(gg)) < 2L) { keep[k] <- FALSE; next }
        #t_j[k] <- compute_AreaPDF(d_loo, target, group, align = "approx")$AreaPDF
        t_j[k] <- metric_fn(d_loo)
      }
    }

    t_j <- t_j[keep]
    if (length(t_j) >= 5L && stats::sd(t_j) > 0) {
      tbar <- mean(t_j)
      num  <- sum((tbar - t_j)^3)
      den  <- 6 * (sum((tbar - t_j)^2)^(3/2))
      a    <- if (den == 0 || !is.finite(den)) 0 else num / den
    } else {
      a <- 0  # fallback if jackknife degenerate
    }

    # Adjusted quantile levels for BCa
    z_alpha <- stats::qnorm(c(alpha, 1 - alpha))
    adj_p <- stats::pnorm(z0 + (z0 + z_alpha) / (1 - a * (z0 + z_alpha)))
    adj_p <- pmin(pmax(adj_p, eps), 1 - eps)  # clamp
    ci_bca <- as.numeric(quantile(bootstrap, adj_p, na.rm = TRUE))

    out$ci_bca <- ci_bca
    out$z0 <- z0
    out$a  <- a
    out$level <- conf
  }

  # ----- NEW: Studentized (bootstrap-t) CI -----
  if (isTRUE(studentized) && is.finite(est)) {
    cat(sprintf("studentized is true"))
    # jackknife SE on original
    jackknife_se <- function(data) {
      n <- nrow(data)
      idx <- if (is.finite(jackknife_max) && jackknife_max < n)
        sort(sample.int(n, jackknife_max)) else seq_len(n)
      t_j <- rep(NA_real_, length(idx))
      for (k in seq_along(idx)) {
        d_loo <- data[-idx[k], , drop = FALSE]
        gg <- d_loo[[group]]
        if (length(unique(gg)) < 2L || any(table(gg) == 0)) next
        t_j[k] <- tryCatch(metric_fn(d_loo), error = function(e) NA_real_)
      }
      t_j <- t_j[is.finite(t_j)]
      if (length(t_j) < 5L) return(NA_real_)
      tbar <- mean(t_j)
      sqrt((length(t_j) - 1) / length(t_j) * sum((t_j - tbar)^2))
    }

    se0 <- jackknife_se(df)

    # build t* using the same resamples you already drew:
    # we need SE for each resample; recompute via jackknife on d
    t_star <- rep(NA_real_, B)
    pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
    for (b in seq_len(B)) {
      # regenerate indices to match your design (or keep & store them earlier)
      i <- c(sample(idx0, length(idx0), TRUE), sample(idx1, length(idx1), TRUE))
      d <- df[i, , drop = FALSE]
      gg <- d[[group]]
      if (length(unique(gg)) < 2L || any(table(gg) == 0)) next
      tb  <- tryCatch(metric_fn(d), error = function(e) NA_real_)
      seb <- jackknife_se(d)
      if (!is.finite(tb) || !is.finite(seb) || seb == 0) next
      t_star[b] <- (tb - est) / seb

      Sys.sleep(0.01)  # simulate work

      utils::setTxtProgressBar(pb, b)
    }
    t_star <- t_star[is.finite(t_star)]
    qL <- stats::quantile(t_star, 1 - alpha, na.rm = TRUE)
    qU <- stats::quantile(t_star,     alpha, na.rm = TRUE)
    out$ci_t <- c(est - qL * se0, est - qU * se0)
    out$se0  <- se0
    out$t_star <- t_star
  }

  out
}

#' Bootstrap Confidence Interval for Relative Fairness Improvement
#'
#' Computes bootstrap confidence intervals for the relative improvement
#' of one prediction method over another using fairness metrics.
#'
#' @param df A data frame.
#' @param pred_baseline Column name of the baseline predictions.
#' @param pred_method Column name of the comparison method.
#' @param group Name of the group column.
#' @param cond Optional conditioning variable.
#' @param metric Fairness metric.
#' @param B Number of bootstrap resamples.
#' @param conf Confidence level.
#' @param seed Random seed.
#' @param compute_bca Logical; compute BCa confidence interval.
#' @param n_grid Number of evaluation grid points.
#' @param align_pdf Alignment method passed to `compute_AreaPDF()`.
#' @param jackknife_max Maximum number of jackknife samples.
#'
#' @return A list containing the estimated improvement and confidence
#' intervals.
#'
#' @export
bootstrap_metric_diff_ci <- function(df, pred_baseline, pred_method, group, cond = NULL,
                         metric = c("AreaCDF","AreaPDF", "AreaPDF_cond"),
                         B = 2000, conf = 0.95, seed = 1, compute_bca = TRUE,
                         n_grid = 1000, align_pdf = "approx",jackknife_max = Inf){

  set.seed(seed)
  metric <- match.arg(metric)
  eps <- 1e-8

  # 1) Select the metric function once (expects data + score column) ---
  metric_fn <- switch(
    metric,
    AreaCDF      = function(data, score_col)
      compute_AreaCDF(data, score_col, group_col=group)$AreaCDF,
    AreaPDF      = function(data, score_col)
      compute_AreaPDF(data, score_col, group_col=group, align = align_pdf)$AreaPDF,
    AreaPDF_cond = function(data, score_col)
      compute_AreaPDF_cond(data, score_col, condition_col = cond, group_col = group)$mean_AreaPDF
  )

  # 2) function to computer the statistic,i.e., metric difference for each bootstrap sample
  stat_pct <- function(d, print=FALSE) {
    b <- metric_fn(d, pred_baseline)
    m <- metric_fn(d, pred_method)
    if (print){
      cat(sprintf("Metric baseline: %.4f; method = %.4f\n", b,m))
    }

    (b - m) / max(b, eps)
  }

  safe_stat_pct <- function(d) {
    # your existing stat; add guards against tiny denominators / non-finite
    val <- stat_pct(d)
    if (!is.finite(val)) return(NA_real_)
    # optional: hard cap to kill numerical explosions
    if (abs(val) > 1e3) return(NA_real_)  # pre-specify this
    val
  }


  # 3) computer the statistics on original sample
  point_pct <- stat_pct(df)
  cat(sprintf("Difference is %.4f\n", point_pct))

  # 4) Stratified indices (assumes binary 0/1 group labels)
  g <- df[[group]]; i0 <- which(g == 0); i1 <- which(g == 1)


  bootstrap <- numeric(B)
  for (b in seq_len(B)) {
    i <- c(sample(i0, length(i0), TRUE),
           sample(i1, length(i1), TRUE))
    d_b <- df[i, , drop = FALSE]
    bootstrap[b] <- tryCatch(safe_stat_pct(d_b), error = function(e) NA_real_)
  }

  alpha <- (1 - conf)/2
  ci_pct_perc <- as.numeric(quantile(bootstrap, c(alpha, 1 - alpha), na.rm = TRUE))
  se_pct <- stats::sd(bootstrap, na.rm = TRUE)

  out <- list(point = point_pct, ci_perc = ci_pct_perc, se = se_pct, bootstrap = bootstrap)
  cat("95% CI (percentile): [", round(ci_pct_perc[1], 3), ", ", round(ci_pct_perc[2], 3), "]\n")
  check_bootstrap_skew(bootstrap)

  # --- optional: BCa CI ---
  if (compute_bca) {
    # BCa on % scale
    R <- length(bootstrap); epsR <- 1/(R + 1)
    p0 <- (sum(bootstrap < point_pct) + 0.5 * sum(bootstrap == point_pct)) / R
    p0 <- min(max(p0, epsR), 1 - epsR)
    z0 <- stats::qnorm(p0)

    # jackknife acceleration for % statistic
    n <- nrow(df)
    jk_idx <- if (is.finite(jackknife_max) && jackknife_max < n) sort(sample.int(n, jackknife_max)) else seq_len(n)
    t_j <- numeric(0)
    for (ii in jk_idx) {
      d_loo <- df[-ii, , drop = FALSE]
      gg <- d_loo[[group]]
      if (length(unique(gg)) < 2L || any(table(gg) == 0)) next
      t_j <- c(t_j, stat_pct(d_loo))
    }
    a <- 0
    if (length(t_j) >= 5L) {
      tbar <- mean(t_j)
      num <- sum((tbar - t_j)^3)
      den <- 6 * (sum((tbar - t_j)^2)^(3/2))
      a <- if (den == 0 || !is.finite(den)) 0 else num / den
    }

    z <- stats::qnorm(c(alpha, 1 - alpha))
    adj_p <- stats::pnorm(z0 + (z0 + z) / (1 - a * (z0 + z)))
    adj_p <- pmin(pmax(adj_p, epsR), 1 - epsR)
    ci_pct_bca <- as.numeric(quantile(bootstrap, adj_p, na.rm = TRUE))

    # one-sided improvement (>% 0)
    zL <- stats::qnorm(alpha)
    adj_pL <- stats::pnorm(z0 + (z0 + zL) / (1 - a * (z0 + zL)))
    adj_pL <- min(max(adj_pL, epsR), 1 - epsR)
    L_pct_bca <- as.numeric(quantile(bootstrap, adj_pL, na.rm = TRUE))
    sig_one_sided <- (L_pct_bca > 0)

    out$z0 <- z0
    out$a <- a
    out$level <- conf
    out$ci_pct_bca <- ci_pct_bca
    out$one_sided_lower_pct_bca <- L_pct_bca
    out$sig_one_sided <- sig_one_sided
    out$bootstrap <- bootstrap
  }
  cat("95% CI (BCa): [", round(ci_pct_bca[1], 3), ", ", round(ci_pct_bca[2], 3), "]\n")

  out
}

#' Bootstrap Distribution Diagnostics
#'
#' Computes summary statistics describing the bootstrap distribution,
#' including skewness and Bowley skewness.
#'
#' @param boot Numeric vector of bootstrap estimates.
#' @param main Plot title (currently unused).
#'
#' @return An invisible list containing the mean, median, standard
#' deviation, Fisher skewness, and Bowley skewness.
#'
#' @export
check_bootstrap_skew <- function(boot, main = "Bootstrap Diagnostic") {

  boot <- as.numeric(boot)
  boot <- boot[is.finite(boot)]

  # 1️⃣ Compute basic statistics
  mean_boot <- mean(boot)
  median_boot <- median(boot)
  sd_boot <- stats::sd(boot)
  skew_boot <- e1071::skewness(boot, type = 2)  # Fisher–Pearson
  q <- quantile(boot, c(0.25, 0.5, 0.75))
  bowley <- as.numeric((q[3] + q[1] - 2*q[2]) / (q[3] - q[1]))

  # 2️⃣ Print numeric summary
  #cat("\n====== Bootstrap Skewness Diagnostics ======\n")
  #cat(sprintf("Mean:    %.4f\n", mean_boot))
  #cat(sprintf("Median:  %.4f\n", median_boot))
  #cat(sprintf("SD:      %.4f\n", sd_boot))
  #cat(sprintf("Skewness (Fisher): %.4f\n", skew_boot))
  #cat(sprintf("Bowley skewness:   %.4f\n", bowley))
  #cat("Interpretation:\n")
  #if (abs(skew_boot) < 0.5) {
  #  cat(" → Roughly symmetric.\n")
  #} else if (abs(skew_boot) < 1) {
  #  cat(" → Moderately skewed.\n")
  #} else {
  #  cat(" → Strongly skewed — BCa interval recommended.\n")
  #}

  # 3️⃣ Visualization
  #par(mfrow = c(1, 2))
  #hist(boot, breaks = 30, probability = TRUE, col = "lightgray",
  #     main = paste(main, "- Histogram"), xlab = "Bootstrap statistic")
  #lines(density(boot), col = "blue", lwd = 2)
  #abline(v = mean_boot, col = "red", lwd = 2)
  #abline(v = median_boot, col = "darkgreen", lwd = 2)
  #legend("topright", legend = c("mean", "median"),
   #      col = c("red", "darkgreen"), lwd = 2, cex = 0.8)

  #qqnorm(boot, main = paste(main, "- Q-Q Plot"))
  #qqline(boot, col = "red", lwd = 2)

  invisible(list(
    mean = mean_boot,
    median = median_boot,
    sd = sd_boot,
    skewness = skew_boot,
    bowley = bowley
  ))
}

