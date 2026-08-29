# Homepage copy — draft for review

Edit freely and send back. Every heading below is a **slot** on the page; the text
under it is what currently renders. Delete, rewrite or leave notes inline — I'll
reconcile whatever comes back.

Where I've flagged **[open]**, I need a decision from you.

---

## 1. Header

**Site mark (top left, small):**

Sam Buckberry

**Nav links:**

Questions · Now · Method · Work · Contact

> **[open]** The mark is the only place your name appears above the fold. Keep,
> change, or remove?

---

## 2. Hero

**H1:**

Epigenetics of Health and Disease

**Statement (large serif, sits directly under the H1):**

Genetic variation, environment and life history converge on the epigenome. We
measure that convergence at population scale, to better understand, predict and
prevent disease.

**Lead paragraph:**

The group works on whole-genome DNA methylation in longitudinal human cohorts —
thousands of people, measured over years — to find where the epigenome registers
cardiometabolic disease before it declares itself, and to build the computational
methods that make that record interpretable.

**Affiliation line (small, muted, below a rule):**

Head of Epigenetics, Black Ochre Data Labs — The Kids Research Institute Australia
and the Australian National University. NHMRC Emerging Leadership Fellow. Adjunct
Senior Lecturer, UWA.

> Links in that line: Black Ochre Data Labs, The Kids Research Institute Australia,
> UWA research repository. ANU is currently unlinked — give me a URL if you want it
> linked. The WA Department of Health HREC role has been dropped from the hero and
> lives on the CV; say if you want it back.

---

## 3. Hero figure

**Caption:**

Illustrative. DNA methylation across a promoter CpG island in two groups from one
cohort. Across the flanking sequence they are indistinguishable; at the island they
separate. The finding is the gap.

**Labels inside the figure:**

- Axis: `mCG / CG`, ticks at 0 / 0.5 / 1.0
- Legend: `unaffected`, `affected`
- Annotation: `Δ 34% methylation`
- Region label: `CpG island`
- Track labels: `transcript`, `8 kb`

> **[open]** What should this figure actually show? Options:
>
> 1. **Promoter CpG island difference** (current) — says "we find differences".
> 2. **mQTL** — methylation stratified by genotype at a variant. Says "we connect
>    genetic variation to the epigenome", which matches research question 1.
> 3. **EWAS Manhattan** — says "we work at genome scale".
> 4. **Methylation drift over follow-up** — says "longitudinal cohort, prediction".
>
> Also: swap in real data, or keep it explicitly illustrative?

---

## 4. Research questions

**Section eyebrow:** Open

**Section heading:** Research questions

**Section note:** What the current programme is pointed at.

---

**Q1 — How do genetic variants shape DNA methylation in health and disease?**

Much of the heritability of common disease sits outside coding sequence, in
regulatory variation we can detect but cannot yet interpret. Mapping how variants
act on the methylome is one route from an association to a mechanism.

**Q2 — How do environment and life history become embedded in the epigenome?**

Exposure, development and circumstance leave marks that genotype alone does not
explain. Separating those contributions — and establishing which are stable,
reversible or inherited — is the harder half of the problem.

> Q1 and Q2 are the split of your original two, which overlapped. Q1 is the
> genetic-control question; Q2 is the exposure question. Merge them back into one
> if you'd rather.

**Q3 — Can epigenetic signatures become predictive and diagnostic tools that change care?**

A biomarker that replicates in a cohort is still a long way from something a
clinician can order. We are interested in the tools that reach individuals,
families and communities, and in what it takes for a health system to act on them.

**Q4 — Can we build sovereign data systems for population-scale genomics?**

Analysis at this scale usually means moving data to compute. We are working on the
alternative: systems where communities retain control and participant privacy is
preserved by the architecture rather than by an agreement about it.

**Q5 — How should research ethics work when data outlives the study?**

Genomic data persists, is reanalysed, and now trains models nobody could describe
at the time of consent. Frameworks designed for single studies are governing
open-ended ones, and the gap is widening.

> **[open]** I drafted Q5 from scratch — you gave me the topic only. It's aimed at
> consent outliving the study, which ties to the GigaScience consent review and
> your HREC service. Re-aim it if your actual interest is Indigenous governance
> specifically, or AI in ethics review, or something else.

**Q6 — How do we engineer AI agents to drive the next generation of genomics analysis?**

Pipelines execute decisions somebody already made. The parts of genomics that
still take judgement are the parts worth handing to an agent — and the parts where
we most need it to show its working.

---

## 5. In progress

**Section eyebrow:** In progress

**Section heading:** What's on the bench

**Section note:** Current, not complete.

---

**PROPHECY methylomes**

Whole-genome DNA methylation across a longitudinal cohort of about 1,400 people,
built toward the first population-scale methylome reference for Indigenous
Australians.

*Meta line:* NHMRC Investigator Grant · 2025–2029

> **[open]** "first population-scale methylome reference for Indigenous Australians"
> is my phrasing, not yours — check it's accurate and something you want claimed.

---

**Biomarkers, into clinic**

Taking methylation signatures for cardiometabolic risk from something that
replicates to something that can be measured, priced and acted on where the need
is greatest.

*Meta line:* MRFF Genomics Health Futures Mission

---

**Sovereignty by construction**

Computational systems where communities keep control of their genomic data because
of how the architecture is built, rather than because of what an agreement says.

*Meta line:* With ALIGN and the National Centre for Indigenous Genomics

> **[open]** Check the collaborator line — I inferred it from the team page.

---

**Agent-assisted analysis**

Building and stress-testing agentic workflows for the interpretive passes of
methylome analysis — and working out how to tell when they are wrong.

*Meta line:* Open, ongoing

> **[open]** Agents now appear three times: here, in Q6, and in the Method section
> below. My recommendation is to cut this item and let Q6 and Method carry it.

---

## 6. Method

**Section eyebrow:** Method

**Section heading:** Analysis you can audit

**Section note:** How the computational half actually runs.

---

Most of this work is code. A population-scale methylome study is a long chain of
judgement calls — coverage thresholds, batch structure, cell-type composition,
which of the ten thousand things that correlate with age you are willing to adjust
for — and every link is a place the answer can quietly go wrong without anything
failing.

So we build the chain to be inspected. Pipelines that record why they did what they
did. Agent-driven workflows that take on the interpretive passes an analyst would
otherwise do by hand: checking an assumption against the data, flagging what
doesn't hold, leaving a trace someone else can follow. Written in the open where
the data governance allows it.

The same constraint applies to governance. Where analysis has to run inside a
sovereign environment, auditability is not a nicety — it is the only way anyone
outside can trust a result they cannot inspect the data behind.

**Pull quote:**

The interesting question isn't whether a model can write the analysis. It's whether
it can tell you the analysis is wrong.

**Diagram labels:**

Cohort → Sequencing → Methylome → Model → Clinic, with a dashed feedback box:
`Agent · checks assumptions, records why`

---

## 7. Selected work

**Section eyebrow:** Selected

**Section heading:** Findings worth the space

**Section note:** Five of about sixty. Chosen for what they showed.

> Each entry leads with what was *found*, not where it was published. Journal names
> are present but visually demoted — no gold highlighting.

---

Reprogrammed stem cells hold a memory of the tissue they came from. Passing them
briefly through a naive state erases it, and the corrected cells behave.

*Buckberry, Liu, Poppe et al. · Nature 620, 863–872 · 2023*

---

Gene regulation in the human prefrontal cortex, from gestation to adulthood,
resolved cell by cell across the longest developmental window in the brain.

*Herring et al. · Cell 185, 4428–4447 · 2022*

---

The methylome a marine sponge builds looks convergently like a vertebrate's —
evidence that neuron-flavoured methylation is a solution evolution has reached more
than once.

*de Mendoza et al. · Nature Ecology & Evolution 3, 1464–1473 · 2019*

---

Circulating epigenomic markers track kidney-disease susceptibility in populations
at high risk of type 2 diabetes — a biomarker measurable from blood, not tissue.

*Khurana et al. · Diabetes Research and Clinical Practice 204, 110918 · 2023*

---

Genomic consent has to be revocable, machine-readable and auditable to be
meaningful at scale. A review of how far the technology can currently take that.

*Oliva et al. · GigaScience 13, giae021 · 2024*

---

**Footer link row:** All publications · CV

> **[open]** Are these the right five? Swap any out. The one-line summaries are my
> reading of each paper — correct anything I've overstated or got wrong.

---

## 8. Contact

**Label:** Get in touch

**Email:** sam.buckberry@thekids.org.au

**Note under the email:**

Happy to hear from people working on methylation at scale, on genomic data
governance, or on agents that do real analytical work.

**Label:** Elsewhere

**Links:** GitHub · Scholar · ORCID · LinkedIn

> **[open]** GitHub is currently pointed at `github.com/SamBuckberry`. Confirm
> that's the profile you want linked — it's the single biggest gap on the current
> site for the positioning you described.

**Portrait:** small, greyscale, in this block rather than the hero.

---

## 9. Colophon (bottom rule)

Design concept — Newsreader & Instrument Sans

Perth & Canberra, Australia

> **[open]** Replace the left-hand line with whatever you want there — or drop the
> row entirely.

---

## Not yet drafted

- **CV page** — structure stays as-is. One change proposed: funding amounts move
  from accent colour to muted grey, so the page stops reading as a ledger.
- **Publications page** — structure stays as-is. Journal names lose the bronze
  highlight for the same reason.

Tell me if you want either rewritten rather than restyled.
