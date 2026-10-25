#!/usr/bin/env Rscript
# make_figS20_one_vs_two_epoch_202610.R
#
# Supplementary Fig S20: constant size (one-epoch) vs two-epoch demography on the
# synonymous SFS of the single-sweep neutral controls (no sweep), folded, n = 100,
# +/- 1 / 3 / 5 Mb. Likelihood-ratio test statistic D = 2 x (LL two - LL one) per
# replicate, with the chi-square (2 df) critical value at alpha = 0.05. The one-epoch
# model is the two-epoch model with nu = 1, so a converged two-epoch fit has
# LL >= one-epoch; rows below that are flagged non-converged and set to D = 0.
# Theme and palettes copied from make_manuscript_figures_202610.R.
#
# Inputs: dadi_analysis_folded/<neutral run>/n100_d<d>mb/{one,two}_epoch/demography/demography_results.csv
#         (fit_demo_one_epoch_neutral_202610.sh, fit_demo_dfe_human_202610.sh)
# Output: manuscript_figures/supplement/FigS20_one_vs_two_epoch.{pdf,png},
#         manuscript_figures/Data_per_replicate_lrt_one_vs_two.csv,
#         manuscript_figures/TableS11_lrt_one_vs_two_epoch.csv
#
# Usage: module load R/4.2.2; Rscript make_figS20_one_vs_two_epoch_202610.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE_DIR <- "/u/home/e/eewade/project-awbigham/data/working-data/202610"
FIT_DIR  <- file.path(BASE_DIR, "dadi_analysis_folded")
OUT_DIR  <- file.path(BASE_DIR, "manuscript_figures")
SUPP_DIR <- file.path(OUT_DIR, "supplement")
N_VAL    <- 100
D_VALS   <- c(1, 3, 5)
D_LEVELS <- paste0("+/- ", D_VALS, " Mb")

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
PAL_R <- c("r = 1e-8" = "#E69F00", "r = 1e-9" = "#0072B2")
rec_lab <- function(r) factor(ifelse(as.numeric(r) > 5e-9, "r = 1e-8", "r = 1e-9"), levels = names(PAL_R))

LRT_DF <- 2
LRT_CRIT <- qchisq(0.95, LRT_DF)
LRT_CRIT_LABEL <- sprintf("Critical value (chi-square, %d df, alpha = 0.05; D = %.2f)", LRT_DF, LRT_CRIT)

runs <- list.dirs(FIT_DIR, recursive = FALSE)
runs <- runs[grepl("^humansim_neutral_dfehuman_", basename(runs))]
ll <- rbindlist(lapply(runs, function(r) {
  m <- regmatches(basename(r), regexec("_rec([^_]+)_seq[0-9]+_rep([0-9]+)$", basename(r)))[[1]]
  rbindlist(lapply(D_VALS, function(d) rbindlist(lapply(c("one_epoch", "two_epoch"), function(mod) {
    f <- file.path(r, sprintf("n%d_d%dmb", N_VAL, d), mod, "demography", "demography_results.csv")
    if (!file.exists(f)) return(NULL)
    data.table(rec = m[2], rep = as.integer(m[3]), d = d, model = mod,
               ll = fread(f)$multinomial_log_likelihood[1])
  }))))
}))
dll <- dcast(ll, rec + rep + d ~ model, value.var = "ll")
n_missing <- dll[!is.finite(one_epoch) | !is.finite(two_epoch), .N]
if (n_missing > 0) warning(n_missing, " replicate x window rows are missing a fit and are dropped")
dll <- dll[is.finite(one_epoch) & is.finite(two_epoch)]
dll[, `:=`(dll = two_epoch - one_epoch, converged = two_epoch - one_epoch >= -1e-6)]
dll[, D := pmax(2 * dll, 0)][, p := pchisq(D, LRT_DF, lower.tail = FALSE)]
dll[, `:=`(rec_f = rec_lab(rec), d_lab = factor(paste0("+/- ", d, " Mb"), levels = D_LEVELS),
           conv = factor(ifelse(converged, "Converged", "Not converged (D = 0)"),
                         levels = c("Converged", "Not converged (D = 0)")))]

fwrite(dll[order(d, rec, rep), .(window = d_lab, rec, rep, ll_one_epoch = one_epoch, ll_two_epoch = two_epoch,
                                 delta_ll = dll, two_epoch_converged = converged, D, p)],
       file.path(OUT_DIR, "Data_per_replicate_lrt_one_vs_two.csv"))
fwrite(dll[, .(replicates = .N, median_delta_ll = signif(median(dll), 3), not_converged = sum(!converged),
               p_lt_0.05 = sum(p < 0.05)), by = .(window = d_lab, rec)][order(window, rec)],
       file.path(OUT_DIR, "TableS11_lrt_one_vs_two_epoch.csv"))

p <- ggplot(dll, aes(x = rec_f, y = D, colour = rec_f)) +
  geom_hline(data = data.frame(y = LRT_CRIT, what = LRT_CRIT_LABEL), aes(yintercept = y, linetype = what),
             colour = "grey40", linewidth = 0.3) +
  scale_linetype_manual(values = setNames("dashed", LRT_CRIT_LABEL), name = NULL) +
  geom_boxplot(aes(fill = rec_f), outlier.shape = NA, width = 0.6, alpha = 0.3, linewidth = 0.3, show.legend = FALSE) +
  geom_jitter(aes(shape = conv), width = 0.15, height = 0, size = 1, alpha = 0.8) +
  facet_grid(. ~ d_lab) +
  scale_colour_manual(values = PAL_R, guide = "none") + scale_fill_manual(values = PAL_R) +
  scale_shape_manual(values = c("Converged" = 16, "Not converged (D = 0)" = 1), name = "Two-epoch fit", drop = FALSE) +
  # Shape legend only when some two-epoch fit did not converge (otherwise every point is filled).
  guides(shape = if (all(dll$converged)) "none" else guide_legend(order = 1), linetype = guide_legend(order = 2)) +
  labs(x = "Recombination rate", y = "LRT statistic D, two-epoch vs one-epoch\n(2 x [LL(two) - LL(one)])") +
  theme_ms + theme(legend.box = "vertical", legend.spacing.y = unit(1, "pt"))

for (ext in c("pdf", "png")) ggsave(file.path(SUPP_DIR, paste0("FigS20_one_vs_two_epoch.", ext)), p,
                                    width = 7.2, height = 3.6, dpi = 300)
message("Saved FigS20_one_vs_two_epoch")
print(fread(file.path(OUT_DIR, "TableS11_lrt_one_vs_two_epoch.csv")))
