# Landing-page figure

`landing-figure.R` generates `assets/figure/landing-figure.svg` — the four-panel
composite on the homepage.

```sh
LC_ALL=C.UTF-8 Rscript figure/landing-figure.R
```

Requires `ggplot2`, `patchwork` and `svglite`.

## The data are simulated

Every panel is generated from `set.seed(19)`, not from real results. The figure
is a schematic of the analysis, and the caption on the page says so. Swapping in
real data means replacing the four data frames — `ewas`, `meth`, `snps` and
`variance` — and leaving the plotting code alone.

## Two things that will bite

**Run it under a UTF-8 locale.** `svglite` silently replaces `−`, `₁₀` and `Δ`
with dots when `LC_CTYPE=C`, producing a figure with mangled axis labels and no
error. The script checks for this and stops rather than writing a broken file.

**Colours are placeholders.** The script draws with literal hex values and then
rewrites them to CSS custom properties (`var(--data-b)` and friends) in the
finished SVG, so one file tracks both the light and dark site themes. If you add
a colour, add it to `token_map` too, or it will be baked in and wrong in one of
the two themes.

## Adjusting the layout

Panel proportions are the `heights` argument to `plot_layout()` at the foot of
the script. The shared locus geometry — window width, where the methylation
change sits, how distal the lead variant is — is the block of constants near the
top (`WINDOW_KB`, `CENTRE_KB`, `ISLAND_HALF`, `MQTL_KB`); panels B, C and the
CG-density track all read from those, so moving the variant moves it everywhere
at once.

Two `theme()` traps worth knowing, both already handled: setting `axis.text` or
`axis.line` in a later `theme()` call does *not* override an `axis.text.y` or
`axis.line.y` set earlier by `theme_locus()`, so those are blanked by name; and
`patchwork` paints its own background, so the composition sets a transparent
`plot.background` with `&` after the layout.

---

# Landing-page schematic

`landing-schematic.py` generates `assets/figure/landing-schematic.svg` — a
Nature Reviews-style concept figure covering the same ground as the data figure,
but explaining the mechanism rather than showing a result.

```sh
python3 figure/landing-schematic.py
```

No dependencies. Styling lives in `assets/figure/schematic.css`, which uses the
site's CSS custom properties, so the schematic tracks both themes without a
second copy.

Three panels:

- **a** the four sources converging on DNA methylation
- **b** a distal regulatory variant setting a CpG island through chromatin looping
- **c** the observed variation attributed between those sources

Written in Python rather than R because it is a drawing, not a plot — there is no
data behind it and nothing to recompute. Every quantity it would otherwise have
to invent has been designed out: panel **c** fans into five equal boxes rather
than a proportioned bar, so the figure makes its point without claiming any
number it cannot support. The only magnitude anywhere is the approximate 6 kb
distance in panel **b**.

Colours come from the same CSS variables as the rest of the site. The accent
(`--data-b`) is reserved for things that carry meaning: methylated CpGs, the
variant, and the chromatin loop. Everything structural is drawn in `--muted` or
`--faint`.
