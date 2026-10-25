#!/usr/bin/env Rscript
# make_manuscript_figures_202610.R
#
# Manuscript figures from the 202610 genic-only resimulation (human DFE of
# Kim et al. 2017 only). Main-text defaults: folded spectra, n = 1000,
# d = 5 Mb (whole 10 Mb chromosome), both r, all sweep s.
#
#   main/Fig1_syn_sfs_distortion     single sweep: A synonymous delta theta_W
#                                    (1 - pi / theta_W); B nu; C syn SFS (observed, constant size,
#                                    two-epoch fit)
#   main/Fig2_dfe_bias               single sweep: A binned DFE (two-epoch); B nonsyn SFS
#                                    (observed, two-epoch + gamma fit, constant size + Kim DFE);
#                                    C -E[s] vs gamma shape, one point per replicate
#   main/Fig3_recurrent_sweeps       recurrent sweeps, every s / beneficial proportion and r:
#                                    A nu; B binned DFE (two-epoch)
#   supplement/FigS1_rare_excess     delta theta_W = 1 - pi / theta_W, syn and nonsyn, every condition
#   supplement/FigS2_syn_sfs         synonymous SFS grid (rows +/- 1/3/5 Mb, columns s x r)
#   supplement/FigS3_syn_sfs_residuals  two-epoch residuals on the synonymous SFS
#   supplement/FigS4_two_vs_three_epoch  nu, LL(three) - LL(two) and LRT statistic D, rows +/- 1/3/5 Mb
#   supplement/FigS5_demography_choices  two-epoch nu by n, window, fold
#   supplement/FigS6_nonsyn_sfs      nonsynonymous SFS grid (same layout as S2)
#   supplement/FigS7_window_distance d = 1 / 3 / 5 Mb (S7, S8, S11: A binned DFE, B -E[s] vs gamma shape)
#   supplement/FigS8_sample_size     n = 100 vs 1000
#   supplement/FigS9, S10            S2 / S6 grids at n = 1000
#   supplement/FigS11_folded_unfolded folded vs unfolded
#   supplement/FigS12, S13           S2 / S6 grids, unfolded
#   supplement/FigS14_three_epoch_dfe three-epoch binned DFEs (single and recurrent)
#   supplement/FigS15_recurrent_sweep_trajectories  fixed beneficials over time
#   supplement/FigS16_recurrent_sfs_dfe_params  recurrent sweeps: syn SFS and -E[s] vs gamma shape
#   supplement/FigS16_recurrent_sfs_rare_excess  alternate S16: syn SFS and delta theta_W (observed vs two-epoch fit)
#   TableS10_snp_counts.csv, TableS3_neutral_syn_sfs_vs_snm.csv, TableS1_simulation_parameters.csv
#
# Inputs:
#   working-data/202610/comparisons/dfe_by_fit.csv
#       (compare_dfe_single_recurrent_202610.py). The bin columns there are in dadi
#       s units; this script re-bins from shape and scale in SLiM units (x 2) so the
#       fits are on the same scale as the Kim et al. 2017 truth.
#   dadi_analysis_folded / dadi_analysis_recurrent_folded  demography fits and spectra
#   manuscript_figures/data/constant_kim_nonsyn_sfs.csv
#       (expected_sfs_constant_kim_202610.py)
#   tracksweep_output / recurrentsweep_output  SLiM summaries and fixed beneficials
#   rare_excess/rare_excess_by_rep.csv  (rare_excess_202610.py; run before this script)
#
# Usage:
#   module load R/4.2.2
#   Rscript make_manuscript_figures_202610.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

BASE_DIR <- "/u/home/e/eewade/project-awbigham/data/working-data/202610"
OUT_DIR  <- file.path(BASE_DIR, "manuscript_figures")
MAIN_DIR <- file.path(OUT_DIR, "main")
SUPP_DIR <- file.path(OUT_DIR, "supplement")
DATA_DIR <- file.path(OUT_DIR, "data")
for (d in c(MAIN_DIR, SUPP_DIR, DATA_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

TRUE_SHAPE     <- 0.186
TRUE_MEAN_S    <- 0.01314833   # SLiM units
TRUTH_LABEL    <- "True human DFE (Kim et al. 2017)"
DADI_TO_SLIM_S <- 2            # dadi fitnesses 1, 1+2hs, 1+2s; SLiM 1, 1+hs, 1+s

MAIN_FOLD <- "folded"
MAIN_N    <- 100
MAIN_D    <- 1   # single-sweep window (+/- Mb around the sweep); recurrent fits are whole-chromosome (5)
REC_D     <- 5

# Lin et al. 2026 (PLoS Biol) Fig 1D bins, in s notation, with the two extremes named.
S_BIN_EDGES  <- c(0, 1e-5, 1e-4, 1e-3, 1e-2, Inf)
S_BIN_LABELS <- c("-1e-5 < s\n<= 0\nNearly\nneutral", "-1e-4 < s\n<= -1e-5", "-1e-3 < s\n<= -1e-4",
                  "-0.01 < s\n<= -1e-3", "s <= -0.01\n\nStrongly\ndeleterious")
discretize_gamma <- function(shape, scale) diff(pgamma(S_BIN_EDGES, shape = shape, scale = scale))
true_bins <- discretize_gamma(TRUE_SHAPE, TRUE_MEAN_S / TRUE_SHAPE)

# ------------------------------------------------------------
# Labels
# ------------------------------------------------------------
S_LEVELS <- c("s = 0 (neutral)", "s = 0.01", "s = 0.05")
REC_LEVELS <- c("r = 1e-8", "r = 1e-9")
# Recurrent conditions: s of each beneficial and the proportion of new nonsyn
# mutations that are beneficial (p * s = 1e-6 in every condition).
REC_COND <- data.table(ben_s = c(0, 0.00032, 0.001, 0.0032, 0.01, 0.05),
                       benfrac = c(NA, 0.003125, 0.001, 0.0003125, 0.0001, 0.00002))
REC_COND[, cond := ifelse(ben_s == 0, "Neutral",
                          sprintf("s = %g\n%g%%\nbeneficial", ben_s, 100 * benfrac))]
COND_LEVELS <- REC_COND$cond

rec_lab <- function(r) factor(ifelse(as.numeric(r) > 5e-9, "r = 1e-8", "r = 1e-9"), levels = REC_LEVELS)
s_lab   <- function(s) factor(ifelse(as.numeric(s) == 0, S_LEVELS[1], paste("s =", as.numeric(s))), levels = S_LEVELS)
cond_lab <- function(s) factor(REC_COND$cond[match(as.numeric(s), REC_COND$ben_s)], levels = COND_LEVELS)
model_lab <- function(m) factor(m, levels = c("two_epoch", "three_epoch"), labels = c("Two-epoch", "Three-epoch"))
fold_lab <- function(f) factor(f, levels = c("folded", "unfolded"), labels = c("Folded", "Unfolded"))
freq_xlab <- function(f) if (f == "folded") "Minor allele frequency" else "Derived allele frequency"
# Sweep s is written s_ben (subscript, plotmath) on every figure: axis titles, and the
# strip, legend and tick labels built from the "s = ..." levels above. Each line of a
# label is a string, except lines starting "s = "; multi-line labels are stacked with atop.
S_BEN_X     <- expression(s[ben])
S_BEN_X_REC <- expression("Recurrent"~s[ben])
s_ben_pm <- function(x, one_line = FALSE) vapply(strsplit(as.character(x), "\n", fixed = TRUE), function(l) {
  if (!length(l)) return('""')
  l <- ifelse(startsWith(l, "s = "), paste0("s[ben] == ", vapply(sub("^s = ", "", l), deparse, "")), vapply(l, deparse, ""))
  if (one_line) return(paste(l, collapse = '*", "*'))   # angled tick labels: atop does not rotate cleanly
  ts <- function(e) paste0("textstyle(", e, ")")
  switch(min(length(l), 3), l, sprintf("atop(%s, %s)", ts(l[1]), ts(l[2])),
         sprintf("atop(%s, textstyle(atop(%s, %s)))", ts(l[1]), ts(l[2]), ts(paste(l[-(1:2)], collapse = "~"))))
}, "")
s_ben_labels <- function(x, one_line = FALSE) parse(text = s_ben_pm(x, one_line))
label_s_ben <- structure(function(labels, multi_line = TRUE) label_parsed(lapply(labels, s_ben_pm), multi_line),
                         class = c("function", "labeller"))
S_BEN_LAB <- labeller(sweep_s = label_s_ben, cond = label_s_ben, sr = label_s_ben)

# ------------------------------------------------------------
# Theme and palettes
# ------------------------------------------------------------
# Same look as the 202609 figures (theme_bw, grey strips, legend at the bottom).
# Figures are saved at print size (7.2 in = full page width), so these point
# sizes are what appears on the page (202609: base 10 on 11-13 in canvases,
# i.e. ~6-7 pt after reduction).
theme_ms <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(size = 8.5, margin = margin(2, 2, 2, 2)),
        axis.title = element_text(size = 9.5),
        axis.text = element_text(size = 7.5),
        legend.position = "bottom", legend.text = element_text(size = 8.5), legend.title = element_text(size = 8.5),
        legend.key.size = unit(10, "pt"), legend.margin = margin(0, 0, 0, 0),
        plot.tag = element_text(face = "bold", size = 12),
        plot.title = element_text(size = 9.5, face = "bold"))
angled_x <- theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7))
dfe_x    <- theme(axis.text.x = element_text(size = 6.5, lineheight = 0.85))

# Colour-blind-safe palettes: Okabe & Ito (2008) for categories, viridis for
# ordered series. Vermillion (OI_VERM) is reserved for truth / reference lines.
OI <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73", yellow = "#F0E442",
        blue = "#0072B2", verm = "#D55E00", purple = "#CC79A7", grey = "#999999")
OI_VERM   <- unname(OI["verm"])
PAL_S     <- c("s = 0 (neutral)" = unname(OI["grey"]), "s = 0.01" = unname(OI["sky"]), "s = 0.05" = unname(OI["blue"]))
PAL_COND  <- setNames(c(unname(OI["grey"]), scales::viridis_pal(begin = 0.15, end = 0.85, direction = -1)(length(COND_LEVELS) - 1)), COND_LEVELS)
PAL_R     <- c("r = 1e-8" = unname(OI["orange"]), "r = 1e-9" = unname(OI["blue"]))
PAL_MODEL <- c("Two-epoch" = unname(OI["green"]), "Three-epoch" = unname(OI["purple"]))
PAL_N     <- c("n = 100" = unname(OI["orange"]), "n = 1000" = unname(OI["blue"]))
PAL_D     <- setNames(scales::viridis_pal(begin = 0.15, end = 0.85, direction = -1)(3), c("1 Mb", "3 Mb", "5 Mb"))
PAL_FOLD  <- c("Folded" = unname(OI["orange"]), "Unfolded" = unname(OI["sky"]))
rec_s_labels <- scale_x_discrete(labels = function(x) ifelse(x == "Neutral", "0", sub("^s = ([^\n]*)\n.*$", "\\1", x)))
PAL_SIM   <- c("Neutral" = unname(OI["grey"]), "Single sweep" = unname(OI["blue"]), "Recurrent sweeps" = unname(OI["purple"]))

save_fig <- function(p, dir, name, width, height) {
  for (ext in c("pdf", "png")) ggsave(file.path(dir, paste0(name, ".", ext)), p,
                                      width = width, height = height, dpi = 300)
  message("Saved ", name)
}
tags <- plot_annotation(tag_levels = "A")

# ------------------------------------------------------------
# Per-fit DFEs (single sweep "202610" and "recurrent")
# ------------------------------------------------------------
CMP <- file.path(BASE_DIR, "comparisons")
dfe <- fread(file.path(CMP, "dfe_by_fit.csv"), colClasses = list(character = c("rec", "ben_s", "benfrac")))
# Recurrent conditions shown in Fig 3 (keyed on s and p; runs added later with
# other (s, p) pairs are left out until Fig 3 is extended to them).
REC_KEYS <- REC_COND[ben_s > 0, paste(ben_s, benfrac)]
in_rec_cond <- function(dt) dt$source != "recurrent" | paste(as.numeric(dt$ben_s), as.numeric(dt$benfrac)) %in% REC_KEYS
dfe <- dfe[source %in% c("202610", "recurrent")]
dfe <- dfe[in_rec_cond(dfe)]
dfe[, scale := DADI_TO_SLIM_S * mean_s / shape]   # SLiM units
dfe[, `:=`(rec_f = rec_lab(rec), model_f = model_lab(model))]
dfe[source == "202610", sweep_s := s_lab(ben_s)]
dfe[, cond := cond_lab(ben_s)]

single <- function(dt, fold_val = MAIN_FOLD, n_val = MAIN_N, d_val = MAIN_D)
  dt[source == "202610" & fold == fold_val & n == n_val & d == d_val]
recurrent <- function(dt, fold_val = MAIN_FOLD, n_val = MAIN_N)
  rbind(single(dt, fold_val, n_val, REC_D)[kind == "neutral"], dt[source == "recurrent" & fold == fold_val & n == n_val],
        fill = TRUE)

cat(sprintf("DFE fits: %d single-sweep, %d recurrent\n", dfe[source == "202610", .N], dfe[source == "recurrent", .N]))

# ------------------------------------------------------------
# Shared plotting pieces
# ------------------------------------------------------------
# Binned DFE: one box per level of `fillvar` in each s bin, truth as red diamonds.
dfe_bins_plot <- function(dt, fillvar, pal, facets = . ~ rec_f) {
  b <- dt[, .(bin = factor(S_BIN_LABELS, levels = S_BIN_LABELS), prop = discretize_gamma(shape, scale)),
          by = eval(unique(c(fillvar, all.vars(facets)[all.vars(facets) != "."], "rep")))]
  truth <- data.table(bin = factor(S_BIN_LABELS, levels = S_BIN_LABELS), prop = true_bins, what = TRUTH_LABEL)
  ggplot(b, aes(x = bin, y = prop, fill = .data[[fillvar]])) +
    geom_boxplot(outlier.shape = NA, width = 0.75, alpha = 0.7, linewidth = 0.25,
                 position = position_dodge(width = 0.85, preserve = "single")) +
    geom_errorbar(data = truth, aes(x = bin, ymin = prop, ymax = prop, linetype = what), inherit.aes = FALSE,
                  width = 0.9, colour = OI_VERM, linewidth = 0.5) +
    facet_grid(facets, labeller = S_BEN_LAB) +
    scale_fill_manual(values = pal, name = NULL, drop = FALSE, labels = s_ben_labels) +
    scale_linetype_manual(values = setNames("dashed", TRUTH_LABEL), name = NULL) +
    labs(x = "Selection coefficient (s)", y = "Proportion of DFE mass") +
    theme_ms + dfe_x
}

# Simulations are at constant size, so the true nu is 1.
NU_TRUE_LABEL <- "Simulated constant size (nu = 1)"
box_by <- function(dt, xvar, yvar, fillvar, pal, xlab, ylab, hline = NULL, hline_label = NULL, log_y = FALSE) {
  p <- ggplot(dt, aes(x = .data[[xvar]], y = .data[[yvar]], fill = .data[[fillvar]]))
  if (!is.null(hline) && is.null(hline_label)) p <- p + geom_hline(yintercept = hline, linetype = "dashed", colour = OI_VERM)
  # labelled reference line: drawn through the linetype scale so it gets a legend key
  if (!is.null(hline_label))
    p <- p + geom_hline(data = data.frame(y = hline, what = hline_label), aes(yintercept = y, linetype = what),
                        colour = OI_VERM) +
      scale_linetype_manual(values = setNames("dashed", hline_label), name = NULL)
  p <- p +
    geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.6, linewidth = 0.3, position = position_dodge(width = 0.8)) +
    geom_point(aes(colour = .data[[fillvar]]), size = 0.6, alpha = 0.6, show.legend = FALSE,
               position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8)) +
    facet_grid(. ~ rec_f) +
    scale_fill_manual(values = pal, name = NULL) +
    scale_colour_manual(values = pal) +
    scale_x_discrete(labels = function(x) sub("s = ", "", sub(" \\(", "\n(", x))) +
    labs(x = xlab, y = ylab) +
    theme_ms
  if (log_y) p <- p + scale_y_log10()
  p
}

# ------------------------------------------------------------
# Demography fits and spectra (folded, n = MAIN_N, d = 5 Mb)
# ------------------------------------------------------------
read_fs_vec <- function(f) {
  x <- readLines(f, warn = FALSE)
  v <- as.numeric(strsplit(trimws(x[2]), "\\s+")[[1]])
  v[as.integer(strsplit(trimws(x[3]), "\\s+")[[1]]) == 1] <- NA
  v
}
snm_shape <- function(n, folded) {
  s <- c(0, 1 / seq_len(n - 1), 0)
  if (folded) { i <- 0:n; s <- ifelse(i < n - i, s + rev(s), ifelse(i == n - i, s, 0)) }
  s
}
# Delta theta_W = 1 - pi / theta_W of an SFS vector (entry i + 1 = count i), over the
# polymorphic entries left unmasked by `mask_from` (the observed SFS, so a fitted
# spectrum uses the same entries). pi = sum i(n - i) eta_i / C(n, 2), theta_W = S / a_n;
# for a folded SFS i runs over minor-allele counts, which gives the same pi. As in
# rare_excess_202610.py. Scale-free, so the fitted SFS needs no rescaling.
RX_KIM_LABEL <- "Two-epoch demography + Kim et al. 2017 DFE"
RX_KINDS <- c("Observed", "Two-epoch fit", RX_KIM_LABEL)
delta_theta_w <- function(v, n, mask_from = v) {
  i <- seq_along(v) - 1
  k <- !is.na(mask_from) & !is.na(v) & i > 0 & i < n
  thw <- sum(v[k]) / sum(1 / seq_len(n - 1))
  1 - (sum(i[k] * (n - i[k]) * v[k]) / (n * (n - 1) / 2)) / thw
}
nums <- function(txt) as.numeric(regmatches(txt, gregexpr("[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?", txt))[[1]])
run_re <- paste0("^(humansim_(sweep|neutral)|recurrentsweep)_dfehuman_mu[^_]+_bens([0-9.]+)",
                 "(?:_benfrac([0-9.]+))?_rec([0-9e.-]+)_seq[0-9]+_rep([0-9]+)$")
const_kim <- {
  f <- file.path(DATA_DIR, "constant_kim_nonsyn_sfs.csv")
  if (file.exists(f)) fread(f) else NULL
}

# One replicate: nu, LL and Na per model; syn SFS (observed, constant size, both
# fits) and nonsyn SFS (observed, two-epoch + gamma fit, constant size + Kim DFE),
# all as proportions of SNPs.
read_run <- function(run_dir, n, d, folded = TRUE) {
  run <- basename(run_dir)
  m <- regmatches(run, regexec(run_re, run, perl = TRUE))[[1]]
  if (length(m) == 0) return(NULL)
  combo <- file.path(run_dir, sprintf("n%d_d%dmb", n, d))
  meta <- data.table(kind = if (m[2] == "recurrentsweep") "recurrent" else m[3], ben_s = as.numeric(m[4]),
                     rec = m[6], rep = as.integer(m[7]))
  emp_f <- file.path(combo, "two_epoch", "demography", "empirical_syn_sfs.fs")
  if (!file.exists(emp_f)) return(NULL)
  emp <- read_fs_vec(emp_f)
  keep <- !is.na(emp) & seq_along(emp) > 1 & seq_along(emp) <= n
  S <- sum(emp[keep])
  snm <- snm_shape(n, folded)
  snm_c <- snm[keep] * S / sum(snm[keep])
  pois <- function(mu, x) sum(-mu + x * log(mu) - lgamma(x + 1))
  fit <- list(); rxl <- list(); sfs <- list(Observed = emp[keep] / S, `Standard neutral model` = snm_c / S)
  for (mod in c("two_epoch", "three_epoch")) {
    rf <- file.path(combo, mod, "demography", "demography_results.csv")
    if (!file.exists(rf)) next
    x <- fread(rf)
    if (nrow(x) == 0) next   # fit still running (file written, no result yet)
    x <- x[1]
    p <- nums(x$best_params)
    fs <- read_fs_vec(file.path(combo, mod, "demography", paste0(mod, "_best_model_sfs_scaled.fs")))
    fit[[mod]] <- data.table(model = mod, nu = if (mod == "two_epoch") p[1] else p[2],
                             ll = x$poisson_log_likelihood, ll_const = pois(snm_c, emp[keep]),
                             Na = x$ancestral_Ne_estimate)
    sfs[[paste(model_lab(mod), "fit")]] <- fs[keep] / sum(fs[keep])
    # Pearson residual of the fit, (data - model) / sqrt(model), on SNP counts
    # (the scaled model SFS is in expected counts), as in 202609.
    if (mod == "two_epoch") {
      sfs$obs_cnt <- emp[keep]; sfs$two_cnt <- fs[keep]
      rxl$syn <- data.table(class = "syn", kind_sfs = c("Observed", "Two-epoch fit"),
                            rare_excess = c(delta_theta_w(emp, n), delta_theta_w(fs, n, emp)))
    }
  }
  # Nonsynonymous SFS
  dfe_dir <- file.path(combo, "two_epoch", "dfe")
  ns <- NULL
  if (file.exists(file.path(dfe_dir, "gamma_dfe_best_model_sfs.fs"))) {
    o <- read_fs_vec(file.path(dfe_dir, "empirical_nonsyn_sfs.fs"))
    e <- read_fs_vec(file.path(dfe_dir, "gamma_dfe_best_model_sfs.fs"))
    k2 <- !is.na(o) & seq_along(o) > 1 & seq_along(o) <= n
    ns <- data.table(allele_count = which(k2) - 1, Observed = o[k2] / sum(o[k2]),
                     `Two-epoch + gamma DFE fit` = e[k2] / sum(e[k2]))
    # nonsyn fit: two-epoch + gamma DFE ("Two-epoch fit" in the delta theta_W panels)
    rxl$nonsyn <- data.table(class = "nonsyn", kind_sfs = c("Observed", "Two-epoch fit"),
                             rare_excess = c(delta_theta_w(o, n), delta_theta_w(e, n, o)))
    # fitted two-epoch demography + Kim et al. 2017 DFE (expected_sfs_fitted_demog_kim_202610.py)
    kim_f <- file.path(dfe_dir, "kim_dfe_fitted_demog_sfs.fs")
    if (file.exists(kim_f))
      rxl$kim <- data.table(class = "nonsyn", kind_sfs = RX_KIM_LABEL,
                            rare_excess = delta_theta_w(read_fs_vec(kim_f), n, o))
    e3_f <- file.path(combo, "three_epoch", "dfe", "gamma_dfe_best_model_sfs.fs")
    if (file.exists(e3_f)) { e3 <- read_fs_vec(e3_f); ns[, `Three-epoch + gamma DFE fit` := e3[k2] / sum(e3[k2])] }
    if (!is.null(const_kim)) {
      n_samp <- n; fold_val <- if (folded) "folded" else "unfolded"
      ck <- const_kim[const_kim$n == n_samp & const_kim$fold == fold_val]   # `n` alone would be the column
      ns[ck, `Constant size + Kim et al. 2017 DFE` := i.prop, on = "allele_count"]
    }
    ns <- cbind(meta, ns)
  }
  list(fit = if (length(fit)) cbind(meta, rbindlist(fit)) else NULL,
       syn = cbind(meta, data.table(allele_count = which(keep) - 1), as.data.table(sfs)),
       nonsyn = ns,
       rx = if (length(rxl)) cbind(meta, rbindlist(rxl)) else NULL)
}

read_runs <- function(root, n, d, folded = TRUE) {
  runs <- list.dirs(root, recursive = FALSE)
  runs <- runs[grepl("dfehuman", basename(runs))]
  # recurrent runs: only the (s, p) pairs in REC_COND
  rec_run <- grepl("^recurrentsweep", basename(runs))
  sp <- sub("^recurrentsweep_dfehuman_mu[^_]+_bens([^_]+)_benfrac([^_]+)_.*$", "\\1 \\2", basename(runs))
  sp <- vapply(strsplit(sp, " "), function(v) if (length(v) == 2) paste(as.numeric(v[1]), as.numeric(v[2])) else "", "")
  runs <- runs[!rec_run | sp %in% REC_KEYS]
  out <- lapply(runs, read_run, n = n, d = d, folded = folded)
  lapply(c(fit = "fit", syn = "syn", nonsyn = "nonsyn", rx = "rx"), function(k) rbindlist(lapply(out, `[[`, k), fill = TRUE))
}
message("Reading demography fits and spectra ...")
sgl <- read_runs(file.path(BASE_DIR, "dadi_analysis_folded"), MAIN_N, MAIN_D)
rcr <- read_runs(file.path(BASE_DIR, "dadi_analysis_recurrent_folded"), MAIN_N, REC_D)
# Neutral controls for the recurrent figures: single-sweep neutral runs at the same
# whole-chromosome window as the recurrent fits.
neut_rec <- read_runs(file.path(BASE_DIR, "dadi_analysis_folded"), MAIN_N, REC_D)
neut_rec <- lapply(neut_rec, function(dt) dt[kind == "neutral"])
label_runs <- function(dt) {
  dt[, rec_f := rec_lab(rec)]
  dt[kind != "recurrent", sweep_s := s_lab(ben_s)]
  dt[, cond := cond_lab(ben_s)]
  if ("model" %in% names(dt)) dt[, model_f := model_lab(model)]
  if ("kind_sfs" %in% names(dt))
    dt[, `:=`(class_f = factor(class, levels = c("syn", "nonsyn"), labels = c("Synonymous", "Nonsynonymous")),
              kind_sfs = factor(kind_sfs, levels = RX_KINDS))]
  dt
}
for (k in names(sgl)) { label_runs(sgl[[k]]); label_runs(rcr[[k]]); label_runs(neut_rec[[k]]) }
# Recurrent figures use the single-sweep neutral replicates as the no-sweep control.
rcr <- lapply(names(rcr), function(k) rbind(neut_rec[[k]], rcr[[k]], fill = TRUE))
names(rcr) <- names(sgl)

# SFS panel: median proportion per allele count, observed IQR as a ribbon.
PAL_SFS <- c("Observed" = "black", "Standard neutral model" = unname(OI["grey"]),
             "Two-epoch fit" = unname(PAL_MODEL["Two-epoch"]), "Three-epoch fit" = unname(PAL_MODEL["Three-epoch"]),
             "Two-epoch + gamma DFE fit" = unname(PAL_MODEL["Two-epoch"]),
             "Three-epoch + gamma DFE fit" = unname(PAL_MODEL["Three-epoch"]),
             "Constant size + Kim et al. 2017 DFE" = unname(OI["orange"]))
LT_SFS <- c("Standard neutral model" = "dashed", "Two-epoch fit" = "solid", "Three-epoch fit" = "solid",
            "Two-epoch + gamma DFE fit" = "solid", "Three-epoch + gamma DFE fit" = "solid",
            "Constant size + Kim et al. 2017 DFE" = "dashed")
# facets: rec_f ~ colvar by default; rowvar replaces rec_f as the row variable
# (rec_f must then be folded into colvar, as in the d x (s, r) grids).
sfs_plot <- function(dt, colvar, ylab, rowvar = "rec_f") {
  long <- melt(dt, id.vars = unique(c("kind", "ben_s", "rec", "rep", "allele_count", "rec_f", rowvar, colvar)),
               measure.vars = intersect(names(PAL_SFS), names(dt)), variable.name = "series", value.name = "prop")
  # For n > 100, counts above 10 are pooled into log-spaced bins (mean proportion
  # per allele count within the bin, per replicate) so sparse high-frequency
  # classes do not show as 0/1-SNP banding. n = 100 is plotted per allele count.
  binned <- max(long$allele_count) > 100
  if (binned) {
    brks <- unique(c(0.5, 1:10 + 0.5, round(10^seq(log10(10.5), log10(max(long$allele_count) + 0.5), length.out = 30)) + 0.5))
    long[, bin := findInterval(allele_count, brks)]
  } else long[, bin := allele_count]
  long <- long[, .(prop = mean(prop), x = exp(mean(log(allele_count)))),
               by = eval(unique(c("series", "bin", "rep", "kind", "ben_s", "rec", "rec_f", rowvar, colvar)))]
  s <- long[!is.na(prop), .(med = median(prop),
                            lo95 = quantile(prop, 0.25), hi95 = quantile(prop, 0.75), x = x[1]),
            by = eval(unique(c("series", "bin", "rec_f", rowvar, colvar)))]
  # Shading: spread of the simulated (observed) SFS across replicates.
  obs <- s[series == "Observed"]
  p <- ggplot(s, aes(x = x, y = med, colour = series)) +
    geom_ribbon(data = obs, aes(ymin = lo95, ymax = hi95, fill = "Simulated, interquartile range"), colour = NA) +
    geom_point(data = obs, size = 0.35) +
    geom_line(data = s[series != "Observed"], aes(linetype = series), linewidth = 0.5) +
    facet_grid(as.formula(paste(rowvar, "~", colvar)), labeller = S_BEN_LAB) +
    scale_colour_manual(values = PAL_SFS, breaks = intersect(names(PAL_SFS), unique(as.character(s$series))),
                        name = NULL) +
    scale_linetype_manual(values = LT_SFS, guide = "none") +
    scale_fill_manual(values = c("Simulated, interquartile range" = "grey80"),
                      name = NULL) +
    scale_y_log10() +
    labs(x = if (binned) "Minor allele count (log scale; binned above 10)" else "Minor allele count", y = ylab) +
    guides(colour = guide_legend(order = 1), fill = guide_legend(order = 2)) +
    theme_ms + theme(legend.box = "vertical")
  if (binned) p <- p + scale_x_log10()
  p
}

# ============================================================
# Main figures (built at print size: 7.2 in wide)
# ============================================================
# nu panels show every fit (no outlier cutoff). Two-epoch nu at n = 100 spans
# 0.80-2.41, so the log axis needs no break.
small_strips <- theme(strip.text.x = element_text(size = 6.5, lineheight = 0.85))

# SLiM per-replicate summary files (header rows repeated by segment are dropped).
read_summary <- function(f) {
  x <- readLines(f, warn = FALSE)
  x <- c(x[1], x[-1][!startsWith(x[-1], "segment")])
  fread(text = paste(x, collapse = "\n"))
}

# Delta theta_W = 1 - pi / theta_W, from the empirical folded
# n = 100 SFS the fits read (rare_excess_202610.py -> rare_excess/rare_excess_by_rep.csv).
# 0 in expectation under the standard neutral model; > 0 when rare variants are
# overrepresented. Single sweeps and neutral controls at +/- MAIN_D Mb, recurrent
# sweeps at the whole chromosome (REC_D).
RX_FILE <- file.path(BASE_DIR, "rare_excess", "rare_excess_by_rep.csv")
if (!file.exists(RX_FILE)) stop("Run rare_excess_202610.py first: ", RX_FILE)
rx <- fread(RX_FILE, colClasses = list(character = c("rec", "benfrac")))
rx <- rx[(kind != "recurrent" & d_mb == MAIN_D) |
         (kind == "recurrent" & d_mb == REC_D & paste(as.numeric(ben_s), as.numeric(benfrac)) %in% REC_KEYS)]
rx[, rec_f := rec_lab(rec)]
rx[, class_f := factor(class, levels = c("syn", "nonsyn"), labels = c("Synonymous", "Nonsynonymous"))]
rx[, sim := factor(fifelse(kind == "recurrent", "Recurrent sweeps", fifelse(kind == "neutral", "Neutral", "Single sweep")),
                   levels = c("Neutral", "Single sweep", "Recurrent sweeps"))]
rx[, cond_h := factor(fifelse(kind == "neutral", "Neutral",
                              fifelse(kind == "sweep", paste0("Single\ns = ", ben_s),
                                      paste0("Recurrent\n", sub("\n.*", "", as.character(cond_lab(ben_s)))))),
                      levels = c("Neutral", "Single\ns = 0.01", "Single\ns = 0.05",
                                 paste0("Recurrent\n", sub("\n.*", "", COND_LEVELS[-1]))))]
RX_YLAB <- "ΔθW (1 − θπ / θW)"
SNM_LABEL <- "Standard neutral model"

# Gamma DFE parameters per replicate: x = scale, y = -E[s] = shape * scale (both SLiM
# units, log axes), one point per replicate coloured by `colvar`; the Kim et al. 2017
# truth as dashed vermillion lines.
LOG_BREAKS <- as.vector(outer(c(1, 3), 10^(-4:2)))
LOG_LABELS <- function(x) sapply(x, format, scientific = FALSE, drop0trailing = TRUE)
es_scale_plot <- function(dt, colvar, pal, facets = . ~ rec_f) {
  truth <- data.table(x = TRUE_MEAN_S / TRUE_SHAPE, y = TRUE_MEAN_S, what = TRUTH_LABEL)
  ggplot(dt, aes(x = scale, y = shape * scale, colour = .data[[colvar]])) +
    geom_vline(data = truth, aes(xintercept = x, linetype = what), colour = OI_VERM, linewidth = 0.4) +
    geom_hline(data = truth, aes(yintercept = y, linetype = what), colour = OI_VERM, linewidth = 0.4) +
    geom_point(size = 1.1, alpha = 0.75) +
    facet_grid(facets, labeller = S_BEN_LAB) +
    scale_x_log10(breaks = LOG_BREAKS, labels = LOG_LABELS) +
    scale_y_log10(breaks = LOG_BREAKS, labels = LOG_LABELS) +
    scale_colour_manual(values = pal, name = NULL, drop = FALSE, labels = s_ben_labels) +
    scale_linetype_manual(values = setNames("dashed", TRUTH_LABEL), name = NULL) +
    guides(colour = guide_legend(order = 1, override.aes = list(size = 2, alpha = 1)), linetype = guide_legend(order = 2)) +
    labs(x = "Gamma scale", y = "\u2212E[s]") +
    theme_ms
}

# Gamma DFE parameters per replicate: x = shape, y = -E[s] (SLiM units, log axis),
# one point per replicate coloured by `colvar`, one panel per sweep s x r; the Kim
# et al. 2017 truth as dashed vermillion crosshairs.
shape_es_plot <- function(dt, colvar, pal, facets = rec_f ~ sweep_s) {
  truth <- data.table(x = TRUE_SHAPE, y = TRUE_MEAN_S, what = TRUTH_LABEL)
  ggplot(dt, aes(x = shape, y = shape * scale, colour = .data[[colvar]])) +
    geom_vline(data = truth, aes(xintercept = x, linetype = what), colour = OI_VERM, linewidth = 0.4) +
    geom_hline(data = truth, aes(yintercept = y, linetype = what), colour = OI_VERM, linewidth = 0.4) +
    geom_point(size = 1.1, alpha = 0.85) +
    facet_grid(facets, labeller = S_BEN_LAB) +
    scale_x_continuous(breaks = c(0.1, 0.3, 0.5)) +
    scale_y_log10(breaks = LOG_BREAKS, labels = LOG_LABELS) +
    scale_colour_manual(values = pal, name = NULL, drop = FALSE, labels = s_ben_labels) +
    scale_linetype_manual(values = setNames("dashed", TRUTH_LABEL), name = NULL) +
    labs(x = "Gamma shape", y = "\u2212E[s]") +
    theme_ms
}

# Delta theta_W of the observed SFS and of the two-epoch fitted SFS (synonymous:
# two-epoch demography; nonsynonymous: two-epoch + gamma DFE), per replicate, from
# read_run. Boxes coloured by `xvar`; observed solid, fit pale with dashed outline.
# Rows synonymous / nonsynonymous; dashed vermillion: standard neutral model (0).
RX_FIT_KINDS <- c("Observed", "Two-epoch fit")
rx_fit_plot <- function(dt, xvar, pal, xlab, facets = class_f ~ rec_f) {
  dt <- dt[is.finite(rare_excess) & kind_sfs %in% RX_FIT_KINDS]
  key_fill <- guide_legend(override.aes = list(fill = "grey45"))
  ggplot(dt, aes(x = .data[[xvar]], y = rare_excess, fill = .data[[xvar]])) +
    geom_hline(data = data.frame(y = 0, what = SNM_LABEL), aes(yintercept = y, linewidth = what),
               linetype = "dashed", colour = OI_VERM) +
    scale_linewidth_manual(values = setNames(0.4, SNM_LABEL), name = NULL) +
    geom_boxplot(aes(alpha = kind_sfs, linetype = kind_sfs), colour = "grey20", outlier.shape = NA, width = 0.75,
                 linewidth = 0.3, position = position_dodge(width = 0.85, preserve = "single")) +
    facet_grid(facets, scales = "free_y") +
    scale_fill_manual(values = pal, guide = "none") +
    scale_alpha_manual(values = setNames(c(0.75, 0.15), RX_FIT_KINDS), breaks = RX_FIT_KINDS, name = NULL) +
    scale_linetype_manual(values = setNames(c("solid", "dashed"), RX_FIT_KINDS), breaks = RX_FIT_KINDS, name = NULL) +
    guides(alpha = key_fill, linetype = key_fill) +
    labs(x = xlab, y = RX_YLAB) +
    theme_ms + theme(legend.key.width = unit(16, "pt"))
}

# ---------------- Fig 1: the sweep distorts the synonymous SFS (two-epoch only)
#   A synonymous delta theta_W = 1 - pi / theta_W (+/- 1 Mb, the
#   SFS window), with the standard-neutral-model expectation (0) as a dashed line;
#   B inferred nu; C synonymous SFS (observed with simulated interquartile range,
#   constant size, two-epoch fit). The binned residuals of the fit are supplement S3.
rx2 <- rx[kind != "recurrent" & class == "syn"]
rx2[, sweep_s := s_lab(ben_s)]
p2rx <- box_by(rx2, "sweep_s", "rare_excess", "sweep_s", PAL_S, S_BEN_X,
               "Synonymous ΔθW\n(1 − θπ / θW)",
               hline = 0, hline_label = SNM_LABEL) +
  guides(fill = "none")
p2a <- box_by(sgl$fit[model == "two_epoch"], "sweep_s", "nu", "model_f", PAL_MODEL, S_BEN_X,
              "Inferred current size (nu)", hline = 1, hline_label = NU_TRUE_LABEL, log_y = TRUE) +
  scale_fill_manual(values = PAL_MODEL, breaks = "Two-epoch", name = NULL) +
  guides(fill = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
  theme(legend.box = "vertical", legend.spacing.y = unit(0, "pt"))
p2b <- sfs_plot(sgl$syn[, !"Three-epoch fit"], "sweep_s", "Proportion of synonymous SNPs")
# S3: Pearson residuals of the two-epoch fit in minor-allele-frequency bins:
# (sum observed - sum model) / sqrt(sum model) per bin and replicate;
# median and interquartile range across replicates. At n = 100 the first two bins are
# singletons and doubletons.
MAF_EDGES <- c(0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5)
MAF_LABELS <- paste0(100 * head(MAF_EDGES, -1), "%-", 100 * MAF_EDGES[-1], "%")
res2 <- sgl$syn[, .(freq_bin = factor(MAF_LABELS[findInterval(allele_count / MAIN_N, MAF_EDGES, left.open = TRUE)],
                                      levels = MAF_LABELS), obs_cnt, two_cnt, sweep_s, rec_f, rep, ben_s, rec)]
res2 <- res2[, .(resid = (sum(obs_cnt) - sum(two_cnt)) / sqrt(sum(two_cnt))), by = .(freq_bin, sweep_s, rec_f, ben_s, rec, rep)]
res2s <- res2[is.finite(resid), .(med = median(resid), lo95 = quantile(resid, 0.25), hi95 = quantile(resid, 0.75)),
              by = .(freq_bin, sweep_s, rec_f)]
p2c <- ggplot(res2s, aes(x = freq_bin, y = med)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  geom_pointrange(aes(ymin = lo95, ymax = hi95), size = 0.15, linewidth = 0.35) +
  facet_grid(rec_f ~ sweep_s, labeller = S_BEN_LAB) +
  labs(x = "Minor allele frequency", y = "Synonymous SFS residual\n(two-epoch fit, Pearson)",
       caption = "Points: median across replicates; bars: interquartile range") +
  theme_ms + angled_x + theme(plot.caption = element_text(size = 7, hjust = 0))
save_fig(p2rx + p2a + p2b + plot_layout(design = "AB\nCC", heights = c(1, 1.15)) + tags,
         MAIN_DIR, "Fig1_syn_sfs_distortion", 7.2, 7.2)
save_fig(p2c, SUPP_DIR, "FigS3_syn_sfs_residuals", 7.2, 4.2)
fwrite(dcast(sgl$fit, ben_s + rec + rep ~ model, value.var = c("nu", "Na", "ll")),
       file.path(OUT_DIR, "Data_per_replicate_nu_ll.csv"))
fwrite(res2s[order(rec_f, sweep_s, freq_bin), .(sweep_s, rec_f, freq_bin, median_resid = signif(med, 3), q25 = signif(lo95, 3), q75 = signif(hi95, 3))],
       file.path(OUT_DIR, "TableS5_syn_residuals.csv"))

# ---------------- Fig 2: the DFE is biased
#   A binned DFE (two-epoch); B nonsynonymous SFS (observed with simulated
#   intervals, two-epoch + gamma fit, constant size + Kim DFE); C -E[s] against
#   gamma shape, one point per replicate (two-epoch)
p3a <- dfe_bins_plot(single(dfe)[model == "two_epoch"], "sweep_s", PAL_S)
p3b <- sfs_plot(sgl$nonsyn[, !"Three-epoch + gamma DFE fit"], "sweep_s", "Proportion of nonsynonymous SNPs")
p3c <- shape_es_plot(single(dfe)[model == "two_epoch"], "sweep_s", PAL_S) +
  scale_colour_manual(values = PAL_S, guide = "none") +
  theme(strip.text.x = element_text(size = 7.5), panel.spacing.x = unit(9, "pt"))
# D: DFE mass in the two extreme bins against window size (+/- 1, 3, 5 Mb around
# the sweep); median and interquartile range across replicates, truth dashed.
dist <- dfe[source == "202610" & fold == MAIN_FOLD & n == MAIN_N & model == "two_epoch",
            { b <- discretize_gamma(shape, scale); .(bin = c("Nearly\nneutral", "Strongly\ndeleterious"), prop = b[c(1, 5)]) },
            by = .(sweep_s, rec_f, d, rep)]
dist <- dist[, .(med = median(prop), lo = quantile(prop, 0.25), hi = quantile(prop, 0.75)), by = .(sweep_s, rec_f, d, bin)]
truth_d <- data.table(bin = c("Nearly\nneutral", "Strongly\ndeleterious"), prop = true_bins[c(1, 5)], what = TRUTH_LABEL)
p3d <- ggplot(dist, aes(x = d, y = med, colour = sweep_s)) +
  geom_hline(data = truth_d, aes(yintercept = prop, linetype = what), colour = OI_VERM, linewidth = 0.4) +
  scale_linetype_manual(values = setNames("dashed", TRUTH_LABEL), name = NULL) +
  geom_line(position = position_dodge(width = 0.4), linewidth = 0.5) +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.4), size = 0.12, linewidth = 0.35) +
  facet_grid(bin ~ rec_f, scales = "free_y") +
  scale_x_continuous(breaks = c(1, 3, 5), labels = c("1", "3", "5")) +
  scale_colour_manual(values = PAL_S, guide = "none") +
  labs(x = "Window around sweep (+/- Mb)", y = "Proportion of DFE mass") +
  theme_ms + theme(legend.key.width = unit(20, "pt"))
# Main Fig 2: A binned DFE, B -E[s] vs gamma shape per replicate, C DFE mass vs window size.
# The nonsynonymous SFS (with constant size + Kim DFE) is supplement S6.
save_fig(p3a + p3c + p3d + plot_layout(design = "AA\nBC", heights = c(1, 1.05), widths = c(1.45, 1)) + tags,
         MAIN_DIR, "Fig2_dfe_bias", 7.2, 7)
# S6 (nonsynonymous) and S2 (synonymous) raw SFS grids are built in the supplement section.
# Alternate Fig 2 keeping the nonsynonymous SFS as panel B (A DFE, B nonsyn SFS,
# C -E[s] vs gamma shape, D window size).
save_fig(p3a + p3b + p3c + p3d + plot_layout(design = "AA\nBB\nCD", heights = c(1, 1.1, 1.05), widths = c(1.45, 1)) + tags,
         MAIN_DIR, "Fig2_dfe_bias_alt_with_nonsyn_sfs", 7.2, 10)

# ---------------- Fig 3: recurrent sweeps
#   A nu; B synonymous SFS; C binned DFE (two-epoch); D -E[s] vs gamma shape per replicate
p4a <- box_by(rcr$fit[model == "two_epoch"], "cond", "nu", "model_f", PAL_MODEL, S_BEN_X_REC,
              "Inferred current size (nu)", hline = 1, hline_label = NU_TRUE_LABEL, log_y = TRUE) + rec_s_labels +
  scale_fill_manual(values = PAL_MODEL, breaks = "Two-epoch", name = NULL)
p4b <- sfs_plot(rcr$syn[, !"Three-epoch fit"], "cond", "Proportion of synonymous SNPs") + small_strips
p4c <- dfe_bins_plot(recurrent(dfe)[model == "two_epoch"], "cond", PAL_COND) +
  guides(fill = guide_legend(nrow = 1, order = 1), linetype = guide_legend(order = 2)) +
  theme(legend.box = "vertical", legend.spacing.y = unit(0, "pt"), legend.text = element_text(size = 7, lineheight = 0.85))
p4d <- shape_es_plot(recurrent(dfe)[model == "two_epoch"], "cond", PAL_COND, facets = . ~ rec_f) +
  scale_x_continuous(breaks = c(0.15, 0.2, 0.25)) +
  guides(colour = guide_legend(nrow = 1, order = 1, override.aes = list(size = 2, alpha = 1))) +
  theme(legend.box = "vertical", legend.spacing.y = unit(0, "pt"), legend.text = element_text(size = 7, lineheight = 0.85))
# Main Fig 3: A nu, B binned DFE. Synonymous SFS and -E[s] vs gamma shape
# for the recurrent runs are supplement S16.
save_fig(p4a + p4c + plot_layout(design = "A\nB", heights = c(0.8, 1)) + tags,
         MAIN_DIR, "Fig3_recurrent_sweeps", 7.2, 5.6)
save_fig(p4b + p4d + plot_layout(design = "A\nB", heights = c(1, 1.1)) + tags,
         SUPP_DIR, "FigS16_recurrent_sfs_dfe_params", 7.2, 7.4)
# Alternate S16: A synonymous SFS; B delta theta_W of the observed and two-epoch fitted SFS.
save_fig(p4b + (rx_fit_plot(rcr$rx, "cond", PAL_COND, S_BEN_X_REC) + rec_s_labels) +
           plot_layout(design = "A\nB", heights = c(1, 1.1)) + tags,
         SUPP_DIR, "FigS16_recurrent_sfs_rare_excess", 7.2, 7.4)
# ============================================================
# Supplement
# ============================================================
# S1: delta theta_W, every condition, synonymous and nonsynonymous
# (data read before the main figures; Fig 1A uses the synonymous single-sweep rows).
# Neutral controls at both windows: +/- MAIN_D Mb (the single-sweep control) and the whole
# chromosome (REC_D, the matched control for the recurrent sweeps), placed beside the recurrent boxes.
# Dashed: standard neutral model (0); dotted grey: +/- MAIN_D Mb neutral-control median per class and r.
rx_s1 <- rbind(rx, fread(RX_FILE, colClasses = list(character = c("rec", "benfrac")))[kind == "neutral" & d_mb == REC_D],
               fill = TRUE)
rx_s1[, `:=`(rec_f = rec_lab(rec),
             class_f = factor(class, levels = c("syn", "nonsyn"), labels = c("Synonymous", "Nonsynonymous")),
             sim = factor(fifelse(kind == "recurrent", "Recurrent sweeps", fifelse(kind == "neutral", "Neutral", "Single sweep")),
                          levels = c("Neutral", "Single sweep", "Recurrent sweeps")))]
NEUT_S1 <- sprintf("Neutral\n+/- %d Mb", c(MAIN_D, REC_D))
rx_s1[, cond_h := factor(fifelse(kind == "neutral", sprintf("Neutral\n+/- %d Mb", as.integer(d_mb)), as.character(cond_h)),
                         levels = c(NEUT_S1[1], levels(rx$cond_h)[2:3], NEUT_S1[2], levels(rx$cond_h)[-(1:3)]))]
neut_rx <- rx[kind == "neutral", .(med = median(rare_excess)), by = .(class_f, rec_f)]
ps1 <- ggplot(rx_s1, aes(x = cond_h, y = rare_excess, fill = sim)) +
  geom_hline(data = data.frame(y = 0, what = SNM_LABEL), aes(yintercept = y, linetype = what), colour = OI_VERM) +
  geom_hline(data = neut_rx, aes(yintercept = med), linetype = "dotted", colour = "grey40") +
  scale_linetype_manual(values = setNames("dashed", SNM_LABEL), name = NULL) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.7) +
  geom_jitter(aes(colour = sim), width = 0.15, size = 0.6, alpha = 0.6, show.legend = FALSE) +
  facet_grid(class_f ~ rec_f) +
  scale_fill_manual(values = PAL_SIM, name = NULL) + scale_colour_manual(values = PAL_SIM) +
  scale_x_discrete(labels = function(x) s_ben_labels(x, one_line = TRUE)) +
  labs(x = NULL, y = RX_YLAB) +
  theme_ms + theme(axis.text.x = element_text(size = 6.5, lineheight = 0.85, angle = 45, hjust = 1))
save_fig(ps1, SUPP_DIR, "FigS1_rare_excess", 7.2, 5.5)
# Table S2 also lists the neutral controls at the whole-chromosome window (REC_D), the
# matched control for the recurrent sweeps (same rows as Fig S1).
rx_tab <- rx_s1[, .(class, kind, ben_s, benfrac, rec, d_mb, rare_excess)]
fwrite(rx_tab[, .(reps = .N, median_rare_excess = signif(median(rare_excess), 3),
                  q25 = signif(quantile(rare_excess, 0.25), 3), q75 = signif(quantile(rare_excess, 0.75), 3)),
              by = .(class, kind, ben_s, benfrac, rec, window_mb = d_mb)][order(class, kind, rec, ben_s, window_mb)],
       file.path(OUT_DIR, "TableS2_rare_excess.csv"))


# S15: recurrent sweep trajectories. SLiM recorded each fixed beneficial's origin
# and fixation tick (fixed_beneficials_*.tsv), not per-generation frequencies, so
# each sweep is drawn as a segment from origin to fixation. A one replicate per
# condition; B fixed beneficials per 10 Mb segment; C sojourn time.
fb_files <- list.files(file.path(BASE_DIR, "recurrentsweep_output"), "^fixed_beneficials_.*dfehuman.*\\.tsv$",
                       full.names = TRUE)
fb <- rbindlist(lapply(fb_files, function(f) {
  run <- sub("^fixed_beneficials_", "", sub("\\.tsv$", "", basename(f)))
  m <- regmatches(run, regexec(run_re, run, perl = TRUE))[[1]]
  x <- readLines(f, warn = FALSE)
  x <- c(x[1], x[-1][!startsWith(x[-1], "segment")])
  d <- if (length(x) > 1) fread(text = paste(x, collapse = "\n")) else
    data.table(segment = integer(), position = integer(), s = numeric(), origin_tick = integer(), fixation_tick = integer())
  if (!paste(as.numeric(m[4]), as.numeric(m[5])) %in% REC_KEYS) return(NULL)
  cbind(data.table(ben_s = as.numeric(m[4]), rec = m[6], rep = as.integer(m[7])), d)
}), fill = TRUE)
fb[, `:=`(rec_f = rec_lab(rec), cond = cond_lab(ben_s), sojourn = fixation_tick - origin_tick)]
segs_per_rep <- rbindlist(lapply(list.files(file.path(BASE_DIR, "recurrentsweep_output"), "^summary_recurrentsweep_.*\\.csv$",
                                            full.names = TRUE), function(f) {
  run <- sub("^summary_", "", sub("\\.csv$", "", basename(f)))
  m <- regmatches(run, regexec(run_re, run, perl = TRUE))[[1]]
  s <- read_summary(f)
  if (!paste(as.numeric(m[4]), as.numeric(m[5])) %in% REC_KEYS) return(NULL)
  data.table(ben_s = as.numeric(m[4]), rec = m[6], rep = as.integer(m[7]), segment = s$segment)
}))
nfix <- fb[, .N, by = .(ben_s, rec, rep, segment)][segs_per_rep, on = .(ben_s, rec, rep, segment)]
nfix[is.na(N), N := 0L]
nfix[, `:=`(rec_f = rec_lab(rec), cond = cond_lab(ben_s))]
BURNIN <- 1e5
ex <- fb[rep == 1]
s5a <- ggplot(ex, aes(y = factor(segment, levels = sort(unique(segment))), yend = factor(segment, levels = sort(unique(segment))))) +
  geom_vline(xintercept = BURNIN, linetype = "dotted", colour = "grey40") +
  geom_segment(aes(x = origin_tick, xend = fixation_tick, colour = cond), linewidth = 1, alpha = 0.8) +
  facet_grid(rec_f ~ cond, labeller = S_BEN_LAB) +
  scale_colour_manual(values = PAL_COND, guide = "none") +
  scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  labs(x = "Generation (dotted: end of burn-in)", y = "Segment (replicate 1)") +
  theme_ms + theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(angle = 40, hjust = 1, size = 6.5),
                   strip.text.x = element_text(size = 6.5))
s5b <- ggplot(nfix, aes(x = cond, y = N, fill = cond)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.7) +
  geom_jitter(width = 0.15, height = 0.1, size = 0.5, alpha = 0.3) +
  facet_grid(. ~ rec_f) + scale_fill_manual(values = PAL_COND, guide = "none") +
  rec_s_labels + labs(x = S_BEN_X_REC, y = "Fixed beneficials per segment") +
  theme_ms + theme(axis.text.x = element_text(size = 6.5, angle = 45, hjust = 1))
s5c <- ggplot(fb, aes(x = cond, y = sojourn, fill = cond)) +
  geom_boxplot(outlier.size = 0.4, alpha = 0.7, width = 0.7) +
  facet_grid(. ~ rec_f) + scale_fill_manual(values = PAL_COND, guide = "none") +
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  rec_s_labels + labs(x = S_BEN_X_REC, y = "Generations from origin\nto fixation") +
  theme_ms + theme(axis.text.x = element_text(size = 6.5, angle = 45, hjust = 1))
save_fig(s5a / (s5b | s5c) + plot_layout(heights = c(1.3, 1)) + tags,
         SUPP_DIR, "FigS15_recurrent_sweep_trajectories", 7.2, 7.5)
fwrite(fb[, .(fixations = .N, median_sojourn = median(sojourn), median_origin = median(origin_tick)),
          by = .(ben_s, rec)][order(rec, ben_s)], file.path(OUT_DIR, "TableS9_recurrent_fixations.csv"))

# S14: three-epoch DFEs (folded, n = MAIN_N, 5 Mb)
s6a <- dfe_bins_plot(single(dfe)[model == "three_epoch"], "sweep_s", PAL_S) + labs(title = "Single sweep")
s6b <- dfe_bins_plot(recurrent(dfe)[model == "three_epoch"], "cond", PAL_COND) + labs(title = "Recurrent sweeps") +
  guides(fill = guide_legend(nrow = 1, order = 1), shape = guide_legend(order = 2)) +
  theme(legend.box = "vertical", legend.text = element_text(size = 7.5, lineheight = 0.85))
save_fig(s6a / s6b + tags, SUPP_DIR, "FigS14_three_epoch_dfe", 7.2, 7.8)

# Demography fits and spectra for the other analysis choices (single sweep).
read_sgl <- function(n, d, folded = TRUE) {
  x <- read_runs(file.path(BASE_DIR, if (folded) "dadi_analysis_folded" else "dadi_analysis"), n, d, folded)
  for (k in names(x)) label_runs(x[[k]])
  x
}
message("Reading single-sweep fits at other n / d / fold ...")
sgl_d <- c(setNames(list(sgl), as.character(MAIN_D)),
           setNames(lapply(setdiff(c(1, 3, 5), MAIN_D), read_sgl, n = MAIN_N), as.character(setdiff(c(1, 3, 5), MAIN_D))))
# Same windows at n = 1000 (folded) and unfolded (n = MAIN_N), for S5 and the S9, S10, S12, S13 grids.
read_sgl_d <- function(n, folded = TRUE) setNames(lapply(c(1, 3, 5), read_sgl, n = n, folded = folded), c("1", "3", "5"))
sgl_d_n1000 <- read_sgl_d(1000)
sgl_d_unf   <- read_sgl_d(MAIN_N, folded = FALSE)
sgl_n1000 <- sgl_d_n1000[[as.character(MAIN_D)]]
sgl_unf   <- sgl_d_unf[[as.character(MAIN_D)]]

# Two-panel supplement: A binned DFE; B -E[s] against gamma shape (as Fig 2B), one point per
# replicate (rows r, columns sweep s), one colour per level of the analysis choice.
# Single sweep, two-epoch.
choice_fig <- function(fdt, var, pal, name) {
  pa <- dfe_bins_plot(fdt, var, pal, facets = rec_f ~ sweep_s)
  pb <- shape_es_plot(fdt, var, pal, facets = rec_f ~ sweep_s) +
    guides(colour = guide_legend(order = 1, override.aes = list(size = 2, alpha = 1)), linetype = guide_legend(order = 2)) +
    theme(axis.text.x = element_text(size = 6.5))
  save_fig(pa / pb + plot_layout(heights = c(1, 0.9)) + tags, SUPP_DIR, name, 7.2, 7.8)
}
# S8: sample size (folded, +/- MAIN_D Mb)
s2f <- dfe[source == "202610" & fold == MAIN_FOLD & d == MAIN_D & model == "two_epoch"]
s2f[, n_lab := factor(paste("n =", n), levels = names(PAL_N))]
choice_fig(s2f, "n_lab", PAL_N, "FigS8_sample_size")
# S7: distance from the sweep (folded, n = MAIN_N)
s3f <- dfe[source == "202610" & fold == MAIN_FOLD & n == MAIN_N & model == "two_epoch"]
s3f[, d_lab := factor(paste(d, "Mb"), levels = names(PAL_D))]
choice_fig(s3f, "d_lab", PAL_D, "FigS7_window_distance")
# S11: folded vs unfolded (n = MAIN_N, +/- MAIN_D Mb)
s4f <- dfe[source == "202610" & n == MAIN_N & d == MAIN_D & model == "two_epoch"]
s4f[, fold_f := fold_lab(fold)]
choice_fig(s4f, "fold_f", PAL_FOLD, "FigS11_folded_unfolded")

# S17-S19: recurrent sweeps under each analysis choice (whole chromosome; the Fig 3
# conditions only, as read_runs and in_rec_cond keep the REC_COND pairs; control =
# single-sweep neutral runs at the same window). One colour per level of `var`.
# A inferred nu; B binned DFE, rows = s bin, x = condition; C -E[s] against gamma shape.
rec_choice_fig <- function(nu_dt, fdt, var, pal, name) {
  pa <- box_by(nu_dt, "cond", "nu", var, pal, S_BEN_X_REC, "Inferred current size (nu)",
               hline = 1, hline_label = NU_TRUE_LABEL, log_y = TRUE) + rec_s_labels
  bins <- fdt[, .(bin = factor(S_BIN_LABELS, levels = S_BIN_LABELS), prop = discretize_gamma(shape, scale)),
              by = c(var, "cond", "rec_f", "rep")]
  pb <- ggplot(bins, aes(x = cond, y = prop, fill = .data[[var]])) +
    geom_hline(data = data.table(bin = factor(S_BIN_LABELS, levels = S_BIN_LABELS), prop = true_bins, what = TRUTH_LABEL),
               aes(yintercept = prop, linetype = what), colour = OI_VERM, linewidth = 0.5) +
    geom_boxplot(outlier.shape = NA, width = 0.75, alpha = 0.7, linewidth = 0.25,
                 position = position_dodge(width = 0.85, preserve = "single")) +
    facet_grid(bin ~ rec_f, scales = "free_y") +
    scale_fill_manual(values = pal, name = NULL) +
    scale_linetype_manual(values = setNames("dashed", TRUTH_LABEL), name = NULL) +
    rec_s_labels +
    labs(x = S_BEN_X_REC, y = "Proportion of DFE mass") +
    theme_ms + theme(strip.text.y = element_text(angle = 0, size = 6.5, lineheight = 0.85))
  pc <- shape_es_plot(fdt, var, pal, facets = rec_f ~ cond) +
    scale_x_continuous(breaks = c(0.15, 0.2, 0.25)) + small_strips +
    guides(colour = guide_legend(order = 1, override.aes = list(size = 2, alpha = 1)), linetype = guide_legend(order = 2)) +
    theme(axis.text.x = element_text(size = 6))
  save_fig(pa / pb / pc + plot_layout(heights = c(0.8, 2, 1)) + tags, SUPP_DIR, name, 7.2, 9.5)
}
# Recurrent runs plus their neutral controls at another fold or n.
read_rec <- function(n, folded = TRUE) {
  r <- read_runs(file.path(BASE_DIR, if (folded) "dadi_analysis_recurrent_folded" else "dadi_analysis_recurrent"), n, REC_D, folded)
  nr <- read_runs(file.path(BASE_DIR, if (folded) "dadi_analysis_folded" else "dadi_analysis"), n, REC_D, folded)
  for (k in names(r)) { label_runs(r[[k]]); label_runs(nr[[k]]) }
  rbind(nr$fit[kind == "neutral"], r$fit, fill = TRUE)
}
message("Reading recurrent fits at other n / fold ...")
rcr_unf_fit   <- read_rec(MAIN_N, folded = FALSE)
rcr_n1000_fit <- read_rec(1000)
# S17: folded vs unfolded (two-epoch, n = MAIN_N)
s17f <- rbind(recurrent(dfe, "folded"), recurrent(dfe, "unfolded"))[model == "two_epoch"]
s17f[, fold_f := fold_lab(fold)]
rec_choice_fig(rbind(rcr$fit[model == "two_epoch"][, fold_f := fold_lab("folded")],
                     rcr_unf_fit[model == "two_epoch"][, fold_f := fold_lab("unfolded")], fill = TRUE),
               s17f, "fold_f", PAL_FOLD, "FigS17_recurrent_folded_unfolded")
# S18: sample size (two-epoch, folded)
s18f <- rbind(recurrent(dfe, MAIN_FOLD, MAIN_N), recurrent(dfe, MAIN_FOLD, 1000))[model == "two_epoch"]
s18f[, n_lab := factor(paste("n =", n), levels = names(PAL_N))]
rec_choice_fig(rbind(rcr$fit[model == "two_epoch"][, n_lab := factor(paste("n =", MAIN_N), levels = names(PAL_N))],
                     rcr_n1000_fit[model == "two_epoch"][, n_lab := factor("n = 1000", levels = names(PAL_N))], fill = TRUE),
               s18f, "n_lab", PAL_N, "FigS18_recurrent_sample_size")
# S19: two- vs three-epoch (folded, n = MAIN_N)
rec_choice_fig(rcr$fit, recurrent(dfe), "model_f", PAL_MODEL, "FigS19_recurrent_two_vs_three_epoch")

# Delta theta_W of the observed and two-epoch fitted SFS, per replicate, for every
# set above (single sweep: n, window, fold; recurrent: whole chromosome).
rx_all <- rbindlist(c(
  lapply(c("1", "3", "5"), function(d) copy(sgl_d[[d]]$rx)[, `:=`(fold = "folded", n = MAIN_N, d = as.numeric(d))]),
  lapply(c("1", "3", "5"), function(d) copy(sgl_d_n1000[[d]]$rx)[, `:=`(fold = "folded", n = 1000, d = as.numeric(d))]),
  lapply(c("1", "3", "5"), function(d) copy(sgl_d_unf[[d]]$rx)[, `:=`(fold = "unfolded", n = MAIN_N, d = as.numeric(d))]),
  list(rcr$rx[kind == "recurrent"][, `:=`(fold = "folded", n = MAIN_N, d = REC_D)])), fill = TRUE)
fwrite(dcast(rx_all, kind + ben_s + rec + rep + fold + n + d + class ~ kind_sfs, value.var = "rare_excess")[
  order(kind, fold, n, d, class, rec, ben_s, rep)], file.path(OUT_DIR, "Data_per_replicate_rare_excess_obs_vs_fit.csv"))

# S6 / S2: raw SFS grids, rows = window around the sweep (+/- 1, 3, 5 Mb), columns =
# sweep s x r (folded, n = MAIN_N). Observed with simulated interquartile range,
# two- and three-epoch fits (synonymous: also constant size; nonsynonymous: also
# constant size + Kim et al. 2017 DFE).
D_LEVELS  <- paste0("+/- ", c(1, 3, 5), " Mb")
SR_LEVELS <- as.vector(outer(c("s = 0", "s = 0.01", "s = 0.05"), REC_LEVELS, paste, sep = "\n"))
grid_dt <- function(k, lst = sgl_d) {
  x <- rbindlist(lapply(names(lst), function(d) lst[[d]][[k]][, d_lab := factor(paste0("+/- ", d, " Mb"), levels = D_LEVELS)]),
                 fill = TRUE)
  x[, sr := factor(paste(sub(" \\(neutral\\)", "", as.character(sweep_s)), rec_f, sep = "\n"), levels = SR_LEVELS)]
}
grid_theme <- list(theme(strip.text.x = element_text(size = 7, lineheight = 0.85), axis.text.x = element_text(size = 6.5)),
                   guides(colour = guide_legend(nrow = 2, order = 1), fill = guide_legend(order = 2)))
save_fig(sfs_plot(grid_dt("nonsyn"), "sr", "Proportion of nonsynonymous SNPs", rowvar = "d_lab") + grid_theme,
         SUPP_DIR, "FigS6_nonsyn_sfs", 7.2, 6)
save_fig(sfs_plot(grid_dt("syn"), "sr", "Proportion of synonymous SNPs", rowvar = "d_lab") + grid_theme,
         SUPP_DIR, "FigS2_syn_sfs", 7.2, 6)
# S9 / S10 / S12 / S13: the same grids at n = 1000 (folded; counts above 10 log-binned) and for
# unfolded spectra (n = MAIN_N; derived allele count).
unf_x <- labs(x = "Derived allele count")
save_fig(sfs_plot(grid_dt("syn", sgl_d_n1000), "sr", "Proportion of synonymous SNPs", rowvar = "d_lab") + grid_theme,
         SUPP_DIR, "FigS9_syn_sfs_n1000", 7.2, 6)
save_fig(sfs_plot(grid_dt("nonsyn", sgl_d_n1000), "sr", "Proportion of nonsynonymous SNPs", rowvar = "d_lab") + grid_theme,
         SUPP_DIR, "FigS10_nonsyn_sfs_n1000", 7.2, 6)
save_fig(sfs_plot(grid_dt("syn", sgl_d_unf), "sr", "Proportion of synonymous SNPs", rowvar = "d_lab") + grid_theme + unf_x,
         SUPP_DIR, "FigS12_syn_sfs_unfolded", 7.2, 6)
save_fig(sfs_plot(grid_dt("nonsyn", sgl_d_unf), "sr", "Proportion of nonsynonymous SNPs", rowvar = "d_lab") + grid_theme + unf_x,
         SUPP_DIR, "FigS13_nonsyn_sfs_unfolded", 7.2, 6)

# S5: demographic inference (two-epoch nu) under each analysis choice, single sweep:
# A sample size (folded, +/- 1 Mb); B window around the sweep (folded, n = MAIN_N);
# C folded vs unfolded (n = MAIN_N, +/- 1 Mb).
nu_choice <- function(lst, var, pal, title) {
  x <- rbindlist(lapply(names(lst), function(l) lst[[l]]$fit[model == "two_epoch"][, (var) := l]))
  x[, (var) := factor(get(var), levels = names(pal))]
  box_by(x, "sweep_s", "nu", var, pal, S_BEN_X, "Inferred current size (nu)",
         hline = 1, hline_label = NU_TRUE_LABEL, log_y = TRUE) + labs(title = title)
}
s10a <- nu_choice(list(`n = 100` = sgl, `n = 1000` = sgl_n1000)[names(PAL_N)], "n_lab", PAL_N, "Sample size")
s10b <- nu_choice(setNames(sgl_d[c("1", "3", "5")], names(PAL_D)), "d_lab", PAL_D, "Window around the sweep (+/- Mb)")
s10c <- nu_choice(list(Folded = sgl, Unfolded = sgl_unf), "fold_f", PAL_FOLD, "Folded vs unfolded")
save_fig(s10a / s10b / s10c + tags, SUPP_DIR, "FigS5_demography_choices", 7.2, 9)
fwrite(rbindlist(list(
  rbindlist(lapply(c("100", "1000"), function(n) (if (n == "100") sgl else sgl_n1000)$fit[, `:=`(n = as.numeric(n), d = MAIN_D, fold = "folded")])),
  rbindlist(lapply(c("3", "5"), function(d) sgl_d[[d]]$fit[, `:=`(n = MAIN_N, d = as.numeric(d), fold = "folded")])),
  sgl_unf$fit[, `:=`(n = MAIN_N, d = MAIN_D, fold = "unfolded")],
  # recurrent sweeps (whole chromosome); their controls are the neutral d = REC_D rows above
  rcr$fit[kind == "recurrent"][, `:=`(n = MAIN_N, d = REC_D, fold = "folded")]), fill = TRUE)[
  , benfrac := vapply(REC_COND$benfrac[match(as.numeric(ben_s), REC_COND$ben_s)], function(b)
      if (is.na(b)) NA_character_ else format(b, scientific = FALSE, drop0trailing = TRUE), "")][
  kind != "recurrent", benfrac := NA][
  , .(fits = .N, median_nu = signif(median(nu), 3), q25 = signif(quantile(nu, 0.25), 3), q75 = signif(quantile(nu, 0.75), 3)),
  by = .(fold, n, d, model, kind, ben_s, benfrac, rec)][order(model, fold, n, d, kind, rec, ben_s)],
  file.path(OUT_DIR, "TableS4_nu_by_choice.csv"))

# S4: two- vs three-epoch demography (single sweep, folded, n = MAIN_N), rows = window
# (+/- 1, 3, 5 Mb). A inferred current size under each model; B LL(three-epoch) -
# LL(two-epoch) on the synonymous SFS per replicate (the models are nested, so a negative
# difference means the three-epoch optimisation did not converge); C likelihood-ratio test
# statistic D = 2 x (LL three - LL two), with the chi-square (2 df) critical value at
# alpha = 0.05 as a dashed line, as in Martinez i Zurita et al. 2025
# (README_lrt_two_vs_three_epoch.md).
LRT_DF <- 2
LRT_CRIT <- qchisq(0.95, LRT_DF)
LRT_CRIT_LABEL <- sprintf("Critical value (chi-square, %d df, alpha = 0.05; D = %.2f)", LRT_DF, LRT_CRIT)
fit_d <- grid_dt("fit")
dll <- dcast(fit_d, kind + ben_s + rec + rep + sweep_s + rec_f + d_lab ~ model, value.var = "ll")
dll[, `:=`(dll = three_epoch - two_epoch, converged = three_epoch - two_epoch >= -1e-6)]
dll[, `:=`(D = pmax(2 * dll, 0))][, p := pchisq(D, LRT_DF, lower.tail = FALSE)]
dll[, conv := factor(ifelse(converged, "Converged", "Not converged (D = 0)"), levels = c("Converged", "Not converged (D = 0)"))]
fwrite(dll[is.finite(dll)][order(d_lab, rec, ben_s, rep),
           .(window = d_lab, ben_s, rec, rep, ll_two_epoch = two_epoch, ll_three_epoch = three_epoch,
             delta_ll = dll, three_epoch_converged = converged, D, p)],
       file.path(OUT_DIR, "Data_per_replicate_lrt.csv"))
fwrite(dll[is.finite(dll), .(replicates = .N, median_delta_ll = signif(median(dll), 3), not_converged = sum(!converged),
                             p_lt_0.05 = sum(p < 0.05)), by = .(window = d_lab, ben_s, rec)][order(window, rec, ben_s)],
       file.path(OUT_DIR, "TableS6_lrt_two_vs_three_epoch.csv"))

no_row_strip <- theme(strip.text.y = element_blank(), strip.background.y = element_blank())
s_x <- scale_x_discrete(labels = function(x) sub(" \\(neutral\\)", "", sub("s = ", "", x)))
s12a <- box_by(fit_d, "sweep_s", "nu", "model_f", PAL_MODEL, S_BEN_X, "Inferred current size (nu)",
               hline = 1, hline_label = NU_TRUE_LABEL, log_y = TRUE) +
  facet_grid(d_lab ~ rec_f) + s_x + no_row_strip +
  guides(fill = guide_legend(order = 1), linetype = guide_legend(order = 2))
s12b <- ggplot(dll[is.finite(dll)], aes(x = sweep_s, y = dll, fill = sweep_s)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.6, linewidth = 0.3) +
  geom_jitter(aes(colour = sweep_s), width = 0.15, height = 0, size = 0.6, alpha = 0.7, show.legend = FALSE) +
  facet_grid(d_lab ~ rec_f) +
  scale_fill_manual(values = PAL_S, guide = "none") + scale_colour_manual(values = PAL_S) + s_x +
  labs(x = S_BEN_X, y = "LL(three-epoch) - LL(two-epoch)\n(synonymous SFS)") +
  theme_ms + no_row_strip
# Non-converged three-epoch fits (LL below two-epoch) have D = 0 (open symbols).
s12c <- ggplot(dll[is.finite(dll)], aes(x = sweep_s, y = D, colour = sweep_s)) +
  geom_hline(data = data.frame(y = LRT_CRIT, what = LRT_CRIT_LABEL), aes(yintercept = y, linetype = what),
             colour = "grey40", linewidth = 0.3) +
  scale_linetype_manual(values = setNames("dashed", LRT_CRIT_LABEL), name = NULL) +
  geom_jitter(aes(shape = conv), width = 0.2, height = 0, size = 0.9, alpha = 0.8) +
  facet_grid(d_lab ~ rec_f) +
  scale_colour_manual(values = PAL_S, guide = "none") +
  scale_shape_manual(values = c(16, 1), name = "Three-epoch fit") + s_x +
  guides(shape = guide_legend(order = 3), linetype = guide_legend(order = 4)) +
  labs(x = S_BEN_X, y = "LRT statistic D, three- vs two-epoch\n(2 x [LL(three) - LL(two)])") +
  theme_ms
save_fig(s12a + s12b + s12c + plot_layout(widths = c(1.15, 1, 1), guides = "collect") + tags &
           theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(1, "pt"),
                 axis.text.x = element_text(size = 6.5)),
         SUPP_DIR, "FigS4_two_vs_three_epoch", 7.2, 7.2)

# Table: SNPs per fit (synonymous and nonsynonymous, unmasked entries of the
# empirical SFS the DFE fit read), median and interquartile range over replicates.
snp_counts <- function(root, n, d, folded) {
  runs <- list.dirs(root, recursive = FALSE)
  runs <- runs[grepl("dfehuman", basename(runs))]
  rbindlist(lapply(runs, function(r) {
    m <- regmatches(basename(r), regexec(run_re, basename(r), perl = TRUE))[[1]]
    if (length(m) == 0) return(NULL)
    if (m[2] == "recurrentsweep" && !paste(as.numeric(m[4]), as.numeric(m[5])) %in% REC_KEYS) return(NULL)
    dd <- file.path(r, sprintf("n%d_d%dmb", n, d), "two_epoch", "dfe")
    if (!file.exists(file.path(dd, "empirical_nonsyn_sfs.fs"))) return(NULL)
    cnt <- function(f) { v <- read_fs_vec(f); sum(v[2:n], na.rm = TRUE) }
    data.table(kind = if (m[2] == "recurrentsweep") "recurrent" else m[3], ben_s = as.numeric(m[4]),
               benfrac = if (nzchar(m[5])) as.numeric(m[5]) else NA_real_, rec = m[6], rep = as.integer(m[7]),
               fold = if (folded) "folded" else "unfolded", n = n, d = d,
               syn = cnt(file.path(dd, "empirical_syn_sfs.fs")), nonsyn = cnt(file.path(dd, "empirical_nonsyn_sfs.fs")))
  }))
}
message("Counting SNPs per fit ...")
snps <- rbindlist(c(
  lapply(seq_len(12), function(i) { g <- expand.grid(n = c(100, 1000), d = c(1, 3, 5), folded = c(TRUE, FALSE))[i, ]
    snp_counts(file.path(BASE_DIR, if (g$folded) "dadi_analysis_folded" else "dadi_analysis"), g$n, g$d, g$folded) }),
  lapply(c(TRUE, FALSE), function(f) lapply(c(100, 1000), function(n)
    snp_counts(file.path(BASE_DIR, if (f) "dadi_analysis_recurrent_folded" else "dadi_analysis_recurrent"), n, REC_D, f))) |> unlist(recursive = FALSE)))
iqr_txt <- function(x) sprintf("%s [%s-%s]", round(median(x)), round(quantile(x, 0.25)), round(quantile(x, 0.75)))
fwrite(snps[, .(reps = .N, syn_snps = iqr_txt(syn), nonsyn_snps = iqr_txt(nonsyn)),
            by = .(kind, ben_s, benfrac, rec, fold, n, d)][order(kind, fold, n, d, rec, ben_s, benfrac)],
       file.path(OUT_DIR, "TableS10_snp_counts.csv"))

# Table: neutral-control synonymous SFS vs the standard neutral model (folded, n = MAIN_N,
# each window). Per replicate, observed / SNM proportion of singletons, of SNPs with
# MAF <= 5% and of SNPs with MAF > 20%; median over replicates. Wilcoxon signed-rank test
# of the singleton ratio against 1 across the 20 independent replicates.
snm_f <- snm_shape(MAIN_N, TRUE)
k_maf <- 2:(MAIN_N / 2 + 1)   # vector entries for minor-allele counts 1..n/2
maf <- (k_maf - 1) / MAIN_N
snm_p <- snm_f[k_maf] / sum(snm_f[k_maf])
neut_snm <- rbindlist(lapply(list.dirs(file.path(BASE_DIR, "dadi_analysis_folded"), recursive = FALSE), function(r) {
  m <- regmatches(basename(r), regexec(run_re, basename(r), perl = TRUE))[[1]]
  if (length(m) == 0 || m[3] != "neutral") return(NULL)
  rbindlist(lapply(c(1, 3, 5), function(d) {
    f <- file.path(r, sprintf("n%d_d%dmb", MAIN_N, d), "two_epoch", "demography", "empirical_syn_sfs.fs")
    if (!file.exists(f)) return(NULL)
    v <- read_fs_vec(f)[k_maf]; obs <- v / sum(v)
    data.table(rec = m[6], rep = as.integer(m[7]), d, S = sum(v),
               single_obs = obs[1], single_ratio = obs[1] / snm_p[1],
               rare_ratio = sum(obs[maf <= 0.05]) / sum(snm_p[maf <= 0.05]),
               common_ratio = sum(obs[maf > 0.2]) / sum(snm_p[maf > 0.2]))
  }))
}))
fwrite(neut_snm[, .(window = paste0("+/- ", d, " Mb"), replicates = .N, median_syn_snps = round(median(S)),
                    singletons_observed = signif(median(single_obs), 3), singletons_snm = signif(snm_p[1], 3),
                    singletons_obs_over_snm = signif(median(single_ratio), 3),
                    replicates_above_snm = sum(single_ratio > 1),
                    wilcoxon_p_singletons = signif(wilcox.test(single_ratio, mu = 1)$p.value, 2),
                    maf_le_5pct_obs_over_snm = signif(median(rare_ratio), 3),
                    maf_gt_20pct_obs_over_snm = signif(median(common_ratio), 3)),
                by = .(d, rec)][order(d, rec)][, d := NULL],
       file.path(OUT_DIR, "TableS3_neutral_syn_sfs_vs_snm.csv"))


# Table: simulation and inference settings (from ConstantPopSizeOneSweep_v4_genic.slim,
# RecurrentSweep_v2_genic.slim, the task-argument files and fit_demo_dfe_human_202610.sh).
fwrite(data.table(rbind(
  c("Population", "Population size (N)", "10,000 diploids, constant", "Wright-Fisher (SLiM)"),
  c("Population", "Burn-in", "100,000 generations (10N)", ""),
  c("Genome", "Chromosome length", "10 Mb", "22 independent segments per replicate; 20 replicates per condition"),
  c("Genome", "Mutation rate (mu)", "1.5e-8 per bp per generation", "exons and introns only"),
  c("Genome", "Recombination rate (r)", "1e-8 or 1e-9 per bp per generation", "uniform"),
  c("Genome", "Intergenic gap", "Uniform(500, 50,000) bp", "no mutations"),
  c("Genome", "Exon length", "lognormal(log 50, log 2) + 1 bp", "30.2% synonymous, 69.8% nonsynonymous"),
  c("Genome", "Intron length", "lognormal(log 100, log 1.5) + 10 bp", "neutral"),
  c("Genome", "Introns per gene", "1 + geometric, continue with probability 0.8", "mean 5 introns, 6 exons"),
  c("Selection", "Nonsynonymous DFE", "gamma, mean s = -0.01314833, shape = 0.186", "Kim et al. 2017"),
  c("Selection", "Dominance (h)", "0.5", "all mutation types"),
  c("Single sweep", "Beneficial s", "0.01, 0.05 (0 = neutral control)", "one mutation at the chromosome centre"),
  c("Single sweep", "Introduction", "generation 100,000, random haplosome", "restarted from the saved state if lost"),
  c("Single sweep", "Sampling time", "at fixation (neutral control: 2,400 generations after introduction)", ""),
  c("Recurrent sweeps", "Beneficial s (p = proportion of nonsynonymous mutations)",
    paste(REC_COND[ben_s > 0, sprintf("%g (p = %g)", ben_s, benfrac)], collapse = "; "), "p x s = 1e-6"),
  c("Recurrent sweeps", "Sampling time", "generation 102,400", ""),
  c("Sampling", "Sample size", "n = 100, 1000 haploid genomes", ""),
  c("Sampling", "SFS windows", "+/- 1, 3, 5 Mb around the centre (single sweep); whole chromosome (recurrent)", ""),
  c("Sampling", "Delta theta_W", "1 - pi / theta_W from the folded n = 100 SFS; +/- 1 Mb (single sweep), whole chromosome (recurrent)", "synonymous and nonsynonymous separately; also for the two-epoch fitted SFS (nonsynonymous: + gamma DFE)"),
  c("Inference", "Demography", "dadi two-epoch and three-epoch on the synonymous SFS", "Poisson likelihood; grid 1200/1600/2000"),
  c("Inference", "DFE", "fitdadi gamma DFE on the nonsynonymous SFS, conditioned on the fitted demography",
    "ns_s_ratio from exact exon/intron bp per window; s converted to SLiM units (x 2)")
)) |> setnames(c("category", "parameter", "value", "notes")), file.path(OUT_DIR, "TableS1_simulation_parameters.csv"))

# ============================================================
# Tables
# ============================================================
bin_cols <- c("nearly_neutral_0_1e-5", "bin2_1e-5_1e-4", "bin3_1e-4_1e-3", "bin4_1e-3_1e-2",
              "strongly_deleterious_gt_1e-2")
bins_median <- function(shape, scale) {
  m <- matrix(t(mapply(discretize_gamma, shape, scale)), ncol = 5)
  as.list(setNames(round(apply(m, 2, median), 3), bin_cols))
}
truth_row <- data.table(source = "truth (Kim et al. 2017)", median_shape = TRUE_SHAPE,
                        t(setNames(round(true_bins, 3), bin_cols)))
tab <- dfe[, c(.(fits = .N, median_shape = signif(median(shape), 3), median_mean_s_slim = signif(median(shape * scale), 3)),
               bins_median(shape, scale)),
           by = .(source, fold, n, d, model, kind, ben_s, benfrac, rec)][order(source, fold, n, d, model, rec, ben_s)]
fwrite(rbind(truth_row, tab, fill = TRUE), file.path(OUT_DIR, "TableS7_dfe_by_condition.csv"))

message("Done. Outputs in ", OUT_DIR)
