# -*- coding: utf-8 -*-
"""Nature Reviews-style concept schematic for the landing page."""
import math

W, H = 640, 458
out = []
A = out.append

def dna(x0, x1, y, rung=11):
    """Two rails with ladder rungs."""
    s = ['<path class="rail" d="M%d %dH%d M%d %dH%d"/>' % (x0, y, x1, x0, y + 7, x1)]
    x = x0 + rung
    while x < x1:
        s.append('<path class="rung" d="M%d %dv7"/>' % (x, y))
        x += rung
    return "".join(s)

def lolli(x, y_dna, filled, h=13, r=3.4):
    cls = "cpg-fill" if filled else "cpg-open"
    return ('<path class="stem" d="M%.1f %dV%.1f"/><circle class="%s" cx="%.1f" cy="%.1f" r="%.1f"/>'
            % (x, y_dna, y_dna - h, cls, x, y_dna - h - r + 1, r))

def box(x, y, w, h, label, cls="unit"):
    return ('<rect class="%s" x="%d" y="%d" width="%d" height="%d" rx="2"/>'
            '<text class="unit-label" x="%.1f" y="%.1f" text-anchor="middle" dy=".33em">%s</text>'
            % (cls, x, y, w, h, x + w / 2, y + h / 2, label))

def title(letter, text, y):
    return ('<text class="pl" x="16" y="%d">%s</text>'
            '<text class="pt" x="30" y="%d">%s</text>' % (y, letter, y, text))

A('<svg class="schem" viewBox="0 0 %d %d" role="img" aria-label="Three-panel schematic. a: genetic variation, development, cell identity and environment converge on DNA methylation. b: a regulatory variant roughly six kilobases away contacts a CpG island through a chromatin loop, so one allele leaves the island unmethylated and the other methylated. c: the observed variation in methylation at that site is attributed between genetics, cell composition, environment, development and an unexplained remainder.">' % (W, H))

# ---------------------------------------------------------------- panel a
A(title("a", "Four sources converge on the epigenome", 16))

inputs = ["Genetic variation", "Development", "Cell identity", "Environment"]
bw, gap = 118, 22
x0 = (W - (4 * bw + 3 * gap)) / 2
centres = []
for i, lab in enumerate(inputs):
    x = x0 + i * (bw + gap)
    centres.append(x + bw / 2)
    A(box(int(x), 28, bw, 24, lab))

MEET = (320, 82)
for cx in centres:
    A('<path class="flow" d="M%.1f 52C%.1f 70 %.1f 64 %.1f %d"/>' % (cx, cx, MEET[0], MEET[0], MEET[1]))
A('<path class="flow" d="M320 82v12"/>')
A('<path class="head" d="M316.5 91l3.5 4 3.5-4"/>')

DNA_A = 102
A(dna(60, 580, DNA_A))
# mixed methylation state across the strand
states = [(84, 1), (116, 0), (148, 1), (180, 1), (212, 0), (250, 1), (286, 0),
          (356, 1), (392, 0), (428, 1), (466, 1), (502, 0), (544, 1)]
for x, f in states:
    A(lolli(x, DNA_A, f))
A('<text class="note" x="60" y="128">DNA methylation</text>')

# ---------------------------------------------------------------- panel b
A(title("b", "A distal variant can set the state of a CpG island", 172))

DNA_B = 240
ISL_X0, ISL_X1 = 150, 250
VAR_X = 470
A('<rect class="band" x="%d" y="%d" width="%d" height="34"/>' % (ISL_X0, DNA_B - 27, ISL_X1 - ISL_X0))
A(dna(60, 580, DNA_B))

A('<text class="note" x="%d" y="%d" text-anchor="middle">CpG island</text>' % ((ISL_X0 + ISL_X1) / 2, DNA_B - 32))

# chromatin loop: physical contact between the two positions
A('<path class="loop" d="M%d %dQ%d 188 %d %d"/>' % (VAR_X, DNA_B - 2, (VAR_X + 200) / 2, 200, DNA_B - 2))
A('<text class="note accent" x="%.1f" y="202" text-anchor="middle">chromatin looping</text>' % ((VAR_X + 200) / 2))
A('<path class="variant" d="M%d %.1fl5 5-5 5-5-5Z"/>' % (VAR_X, DNA_B - 1.5))
A('<text class="note accent" x="%d" y="%d">regulatory variant</text>' % (VAR_X + 12, DNA_B - 10))

# distance between the two
A('<path class="measure" d="M200 %dv6M%d %dv6M200 %dH%d"/>' % (DNA_B + 11, VAR_X, DNA_B + 11, DNA_B + 14, VAR_X))
A('<text class="note" x="%.1f" y="%d" text-anchor="middle">&#8776; 6 kb</text>' % ((200 + VAR_X) / 2, DNA_B + 30))

# the two outcomes
for row, (lab, filled, state) in enumerate(
        [("allele A", False, "island stays unmethylated"),
         ("allele B", True,  "island becomes methylated")]):
    y = 288 + row * 26
    A('<text class="row-label" x="130" y="%d" text-anchor="end" dy=".33em">%s</text>' % (y, lab))
    for i in range(8):
        cls = "cpg-fill" if filled else "cpg-open"
        A('<circle class="%s" cx="%d" cy="%d" r="3.4"/>' % (cls, 146 + i * 13, y))
    A('<text class="note" x="268" y="%d" dy=".33em">%s</text>' % (y, state))

# ---------------------------------------------------------------- panel c
A(title("c", "Attributing the variation", 348))

A('<rect class="observed" x="160" y="364" width="320" height="22" rx="2"/>')
A('<text class="unit-label" x="320" y="375" text-anchor="middle" dy=".33em">variation in methylation at this site</text>')

sources = ["Genetics", "Cell composition", "Environment", "Development", "Unexplained"]
sw, sgap = 104, 20
sx0 = (W - (5 * sw + 4 * sgap)) / 2
for i, lab in enumerate(sources):
    x = sx0 + i * (sw + sgap)
    A('<path class="flow" d="M320 386C320 404 %.1f 398 %.1f 416"/>' % (x + sw / 2, x + sw / 2))
    A('<path class="head" d="M%.1f 412l3.5 4 3.5-4"/>' % (x + sw / 2 - 3.5))
    cls = "unit dashed" if lab == "Unexplained" else "unit"
    A(box(int(x), 420, sw, 24, lab, cls))

A('</svg>')

open("/home/user/sambuckberry.github.io/assets/figure/landing-schematic.svg", "w", encoding="utf-8").write("\n".join(out))
print("wrote schematic:", sum(len(o) for o in out), "chars")
