#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Landing-page figure for sambuckberry.me
#
# Four panels, all illustrative (simulated) data:
#   A  epigenome-wide scan, one locus above threshold
#   B  DNA methylation across that locus in two groups
#   C  association of local variants with methylation at the target CpG
#   D  partition of methylation variance between sources
#
# Between C and D sits a CG-density track carrying the target CpG and a
# sashimi-style arc to the distal lead variant.
#
# Writes a theme-aware SVG: colours are emitted as placeholder hex values and
# then rewritten to CSS custom properties, so one file renders correctly in
# both the light and dark site themes.
#
#   Rscript figure/landing-figure.R
#
# Requires: ggplot2, patchwork, svglite
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(svglite)
})

# svglite silently replaces non-ASCII glyphs (−, ₁₀, Δ) with dots when LC_CTYPE
# is not UTF-8. Fail loudly rather than shipping a figure with mangled labels.
if (!grepl("UTF-8", Sys.getlocale("LC_CTYPE"), fixed = TRUE)) {
  stop("LC_CTYPE is '", Sys.getlocale("LC_CTYPE"), "', not UTF-8.\n",
       "Re-run as:  LC_ALL=C.UTF-8 Rscript figure/landing-figure.R")
}

set.seed(19)

OUT <- "assets/figure/landing-figure.svg"
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)

# --- placeholder colours -----------------------------------------------------
# Each is swapped for a CSS variable after rendering (see token_map below), so
# the exact values here only need to be distinct and roughly theme-correct.
INK   <- "#1B1A17"   # --ink        primary trace, panel letters
MUTED <- "#5B584E"   # --muted      axis titles, secondary points
FAINT <- "#8B8779"   # --faint      tick labels, rules, unexplained segment
RULE  <- "#D8D2C3"   # --rule       axis lines, threshold lines
DATA  <- "#B4562A"   # --data-b     the accent: second group, lead variant, arc
BAND  <- "#F3E4DA"   # --band       shaded CpG-island region
PAPER <- "#F2EFE7"   # --paper      background, segment separators

token_map <- c(
  INK   = "var(--data-a)",
  MUTED = "var(--muted)",
  FAINT = "var(--faint)",
  RULE  = "var(--rule)",
  DATA  = "var(--data-b)",
  BAND  = "var(--band)",
  PAPER = "var(--paper)"
)

# --- shared geometry ---------------------------------------------------------
WINDOW_KB   <- 20          # width of the locus view
CENTRE_KB   <- 10          # position of the methylation change
ISLAND_HALF <- 2           # half-width of the CpG island, kb
MQTL_KB     <- 16          # lead variant: distal, deliberately off-centre
DELTA_LABEL <- "Δ 34%"

base_size <- 9

# A stripped theme: no grid, no panel background. The page supplies the ruled
# paper behind the figure, so drawing our own would double it up.
theme_locus <- function(show_x = FALSE) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid       = element_blank(),
      panel.background = element_blank(),
      plot.background  = element_blank(),
      plot.margin      = margin(2, 4, 2, 4),
      axis.line.y      = element_line(colour = RULE, linewidth = 0.3),
      axis.line.x      = if (show_x) element_line(colour = RULE, linewidth = 0.3) else element_blank(),
      axis.ticks       = element_blank(),
      axis.text.y      = element_text(colour = FAINT, size = base_size - 2),
      axis.text.x      = element_blank(),
      axis.title       = element_blank(),
      plot.title       = element_text(colour = MUTED, size = base_size - 2.4,
                                      face = "plain", hjust = 0,
                                      margin = margin(b = 3)),
      plot.tag         = element_text(colour = INK, size = base_size,
                                      face = "bold", hjust = 0),
      plot.tag.position = c(0, 1)
    )
}

# ============================================================================
# A. Epigenome-wide scan
# ============================================================================
chrom_len  <- c(8.4, 8.2, 6.7, 6.4, 6.1, 5.8, 5.4, 4.9, 4.7, 4.5, 4.6, 4.5,
                3.9, 3.6, 3.4, 3.0, 2.8, 2.7, 2.0, 2.2, 1.6, 1.7)
chrom_end  <- cumsum(chrom_len)
chrom_start <- c(0, head(chrom_end, -1))
chrom_mid  <- (chrom_start + chrom_end) / 2
PEAK_CHR   <- 8

ewas <- do.call(rbind, lapply(seq_along(chrom_len), function(i) {
  n <- max(9, round(chrom_len[i] * 3.4))
  p <- rexp(n, rate = 0.95)
  p[p > 4.6] <- runif(sum(p > 4.6), 1.4, 4.4)
  data.frame(pos = runif(n, chrom_start[i], chrom_end[i]),
             logp = pmin(p, 4.6),
             shade = i %% 2, hit = FALSE)
}))

peak_x <- chrom_start[PEAK_CHR] + chrom_len[PEAK_CHR] * 0.46
peak <- data.frame(
  pos  = peak_x + c(-0.09, -0.055, -0.026, 0, 0.028, 0.059, 0.095, 0.13),
  logp = c(4.9, 6.2, 7.4, 9.3, 7.9, 6.0, 5.1, 4.7),
  shade = 1, hit = TRUE
)
ewas <- rbind(ewas, peak)

bands <- data.frame(xmin = chrom_start, xmax = chrom_end)[seq(2, 22, 2), ]

pA <- ggplot() +
  geom_rect(data = bands,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = BAND, alpha = 0.45) +
  geom_hline(yintercept = 5, colour = FAINT, linetype = "22", linewidth = 0.3) +
  geom_point(data = subset(ewas, !hit),
             aes(pos, logp), colour = MUTED, alpha = 0.5, size = 0.45) +
  geom_point(data = subset(ewas, hit),
             aes(pos, logp), colour = DATA, size = 0.75) +
  annotate("text", x = peak_x, y = 10.3, label = "chr8",
           colour = DATA, size = 2.7, hjust = 0.5) +
  scale_x_continuous(expand = expansion(0, 0)) +
  scale_y_continuous(limits = c(0, 11), expand = expansion(0, 0)) +
  labs(title = "epigenome-wide scan  ·  −log₁₀ P", tag = "A") +
  theme_locus() +
  theme(axis.text.y = element_blank(), axis.line.y = element_blank())

# ============================================================================
# B. Methylation across the locus
# ============================================================================
# Same shape for both groups across the flanks; they part only at the island,
# where one group retains methylation the other does not.
methyl_series <- function(depth, phase) {
  x <- seq(0, WINDOW_KB, length.out = 400)
  z <- (x - CENTRE_KB) / 1.25
  m <- 0.885 - depth * exp(-z^2) +
    0.024 * sin(x * 1.7 + phase) + 0.014 * sin(x * 5.1 + phase * 2.3) +
    runif(length(x), -0.009, 0.009)
  data.frame(x = x, m = pmin(pmax(m, 0.02), 0.97))
}
unaff <- methyl_series(0.795, 0.0)
aff   <- methyl_series(0.455, 1.9)

meth <- rbind(
  transform(unaff, group = "unaffected"),
  transform(aff,   group = "affected")
)
meth$group <- factor(meth$group, levels = c("unaffected", "affected"))

gap <- data.frame(x = unaff$x, lo = unaff$m, hi = aff$m)
gap <- subset(gap, abs(x - CENTRE_KB) <= 2.35)

island <- data.frame(xmin = CENTRE_KB - ISLAND_HALF, xmax = CENTRE_KB + ISLAND_HALF)
y_un <- min(unaff$m); y_af <- min(aff$m)

pB <- ggplot() +
  geom_rect(data = island, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = BAND, alpha = 0.5) +
  geom_hline(yintercept = 0.5, colour = RULE, linetype = "22", linewidth = 0.3) +
  geom_ribbon(data = gap, aes(x = x, ymin = lo, ymax = hi),
              fill = DATA, alpha = 0.16) +
  geom_line(data = meth, aes(x, m, colour = group), linewidth = 0.4) +
  # the one quantity the figure states outright
  annotate("segment", x = CENTRE_KB, xend = CENTRE_KB, y = y_un, yend = y_af,
           colour = DATA, linewidth = 0.4) +
  annotate("segment", x = CENTRE_KB - 0.22, xend = CENTRE_KB + 0.22,
           y = c(y_un, y_af), yend = c(y_un, y_af), colour = DATA, linewidth = 0.3) +
  annotate("segment", x = CENTRE_KB, xend = CENTRE_KB + 3.1,
           y = mean(c(y_un, y_af)), yend = mean(c(y_un, y_af)),
           colour = DATA, linetype = "22", linewidth = 0.3) +
  annotate("text", x = CENTRE_KB + 3.3, y = mean(c(y_un, y_af)), label = DELTA_LABEL,
           colour = DATA, size = 2.7, hjust = 0) +
  scale_colour_manual(values = c(unaffected = INK, affected = DATA), name = NULL) +
  scale_x_continuous(limits = c(0, WINDOW_KB), expand = expansion(0, 0)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                     labels = c("0", "0.5", "1.0"), expand = expansion(0, 0)) +
  labs(title = "mCG / CG", tag = "B") +
  theme_locus() +
  theme(
    legend.position     = c(0.055, 0.30),
    legend.justification = c(0, 0.5),
    legend.key.size     = unit(7, "pt"),
    legend.key          = element_blank(),
    legend.background   = element_blank(),
    legend.text         = element_text(colour = FAINT, size = base_size - 2.4),
    legend.spacing.y    = unit(1, "pt")
  )

# ============================================================================
# C. Variant association at the same locus
# ============================================================================
# The association peak sits several kb from the methylation change: the variant
# acts on the target CpG at a distance.
n_snp <- 58
snp_pos <- runif(n_snp, 0, WINDOW_KB)
d <- abs(snp_pos - MQTL_KB)
snp_logp <- rexp(n_snp, rate = ifelse(d > 2.6, 1.9, 0.75)) *
  (1 + 2.6 * exp(-(d / 1.75)^2))
snps <- data.frame(pos = snp_pos, logp = pmin(snp_logp, 11.4))
lead <- data.frame(pos = MQTL_KB, logp = 11.2)

pC <- ggplot() +
  geom_rect(data = island, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = BAND, alpha = 0.5) +
  geom_hline(yintercept = 5, colour = FAINT, linetype = "22", linewidth = 0.3) +
  geom_point(data = snps, aes(pos, logp), colour = MUTED, alpha = 0.55, size = 0.55) +
  geom_point(data = lead, aes(pos, logp), colour = DATA, size = 1.5, shape = 18) +
  annotate("segment", x = MQTL_KB - 0.25, xend = MQTL_KB - 3.3, y = 11.2, yend = 11.2,
           colour = DATA, linetype = "22", linewidth = 0.3) +
  annotate("text", x = MQTL_KB - 3.5, y = 11.2, label = "lead mQTL",
           colour = DATA, size = 2.7, hjust = 1) +
  scale_x_continuous(limits = c(0, WINDOW_KB), expand = expansion(0, 0)) +
  scale_y_continuous(limits = c(0, 12.6), expand = expansion(0, 0)) +
  labs(title = "variant association  ·  −log₁₀ P", tag = "C") +
  theme_locus() +
  theme(axis.text.y = element_blank(), axis.line.y = element_blank())

# ============================================================================
# CG density track, with the sashimi arc from lead variant to target CpG
# ============================================================================
cg <- sort(c(
  runif(18, 0, CENTRE_KB - ISLAND_HALF),
  runif(22, CENTRE_KB - ISLAND_HALF, CENTRE_KB + ISLAND_HALF),  # island: dense
  runif(16, CENTRE_KB + ISLAND_HALF, WINDOW_KB)
))
cg <- cg[abs(cg - CENTRE_KB) > 0.12]

pTrack <- ggplot() +
  annotate("segment", x = cg, xend = cg, y = 0.72, yend = 1,
           colour = FAINT, linewidth = 0.3) +
  annotate("segment", x = CENTRE_KB, xend = CENTRE_KB, y = 0.72, yend = 1,
           colour = DATA, linewidth = 0.7) +
  # geom_curve gives the sashimi arc in one call; curvature is negative so the
  # arc bows downward, away from the density ticks.
  geom_curve(aes(x = MQTL_KB, xend = CENTRE_KB, y = 0.6, yend = 0.6),
             curvature = -0.30, ncp = 24, colour = DATA, linewidth = 0.4) +
  annotate("text", x = mean(c(MQTL_KB, CENTRE_KB)), y = -0.34,
           label = paste0(MQTL_KB - CENTRE_KB, " kb"),
           colour = DATA, size = 2.7, vjust = 0) +
  annotate("text", x = WINDOW_KB, y = -0.34, label = paste0(WINDOW_KB, " kb"),
           colour = FAINT, size = 2.7, hjust = 1, vjust = 0) +
  scale_x_continuous(limits = c(0, WINDOW_KB), expand = expansion(0, 0)) +
  scale_y_continuous(limits = c(-0.42, 1.05), expand = expansion(0, 0)) +
  labs(title = "CG density") +
  theme_locus() +
  theme(axis.text.y = element_blank(), axis.line.y = element_blank())

# ============================================================================
# D. Where the variance comes from
# ============================================================================
# Proportions are illustrative and deliberately unlabelled: the panel makes a
# qualitative point about attribution, not a quantitative claim.
variance <- data.frame(
  source = factor(
    c("Genetics", "Cell composition", "Environment", "Development", "Unexplained"),
    levels = c("Genetics", "Cell composition", "Environment", "Development", "Unexplained")
  ),
  share = c(34, 21, 18, 12, 15)
)
var_fill <- c(
  Genetics           = DATA,
  `Cell composition` = DATA,
  Environment        = DATA,
  Development        = DATA,
  Unexplained        = FAINT
)
var_alpha <- c(0.92, 0.68, 0.46, 0.28, 0.28)

pD <- ggplot(variance, aes(x = share, y = 1, fill = source, alpha = source)) +
  # orientation = "y" stacks along x, giving one horizontal bar
  geom_col(orientation = "y", width = 0.16, colour = PAPER, linewidth = 0.5,
           position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = var_fill, name = NULL) +
  scale_alpha_manual(values = setNames(var_alpha, levels(variance$source)),
                     name = NULL) +
  scale_x_continuous(expand = expansion(0, 0)) +
  scale_y_continuous(limits = c(0.8, 1.2), expand = expansion(0, 0)) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(colour = NA)),
         alpha = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(title = "variance in methylation explained", tag = "D") +
  theme_locus() +
  theme(
    axis.text.x      = element_blank(),
    axis.text.y      = element_blank(),
    axis.line.x      = element_blank(),
    axis.line.y      = element_blank(),
    legend.position  = "bottom",
    legend.key.size  = unit(6, "pt"),
    legend.text      = element_text(colour = FAINT, size = base_size - 2.4),
    legend.margin    = margin(t = 2),
    legend.box.spacing = unit(2, "pt")
  )

# ============================================================================
# Compose and write
# ============================================================================
fig <- pA / pB / pC / pTrack / pD +
  plot_layout(heights = c(1.05, 1.5, 0.85, 0.7, 0.7)) &
  theme(plot.background = element_rect(fill = "transparent", colour = NA))

tmp <- tempfile(fileext = ".svg")
svglite(tmp, width = 6.6, height = 5.9, bg = "transparent")
print(fig)
invisible(dev.off())

svg <- readLines(tmp, warn = FALSE)

# Swap the placeholder hex values for CSS custom properties so a single file
# tracks the site's light and dark themes. svglite writes colours as lowercase
# hex, sometimes with an alpha suffix, so match case-insensitively and keep any
# trailing alpha pair off the substitution.
for (nm in names(token_map)) {
  hex <- get(nm)
  svg <- gsub(paste0("(?i)", hex), token_map[[nm]], svg, perl = TRUE)
}

writeLines(svg, OUT)
cat("wrote ", OUT, " (", length(svg), " lines)\n", sep = "")
