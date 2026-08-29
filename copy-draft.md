# Homepage copy — round 2

Your edits are applied and live at the prototype. This file is the current state
of the copy. Same convention: each heading is a **slot** on the page, the text
under it is what renders, and **[open]** marks something I still need from you.

Resolved from your last pass is listed at the bottom so you can check I didn't
miss anything.

---

## 1. Header

**Site mark (top left, small):** Sam Buckberry

**Nav links:** Questions · Now · Method · Work · Contact

---

## 2. Hero

**H1:**

Epigenetics of Health and Disease

**Statement (large serif, under the H1):**

Genetic variation, development, cell identity, environment and life history are
all written into the epigenome. We work to resolve their separate contributions,
and to use what that resolves into — to understand, predict and prevent disease.

> **[open]** You asked for alternatives to "converge/convergence" and said the
> measuring didn't quite make sense. The line above is my current pick, but the
> tail ("what that resolves into") is still awkward. Four alternatives:
>
> **A.** The epigenome records the combined action of genetic variation,
> development, cell identity and environment. We work to separate those
> contributions, and to use what they reveal to understand, predict and prevent
> disease.
>
> **B.** Genetic variation, development, cell identity and environment each leave
> their trace on the epigenome. Disentangling those traces is how we get from an
> observed difference to an explanation, and from an explanation to a prediction.
>
> **C.** The epigenome integrates genetic variation, development, cell identity
> and environment. We work to decompose that integration across scales, from
> single cells to whole cohorts, to understand, predict and prevent disease.
>
> **D.** Genetic variation, development, cell identity, environment and life
> history are all inscribed in the epigenome. Our work is to tell them apart —
> and to make what they encode clinically useful.
>
> My preference is **A**: "records the combined action" avoids the convergence
> problem, and "separate those contributions" is what you actually do. **B** is
> the most distinctive if you want something less standard.

**Lead paragraph:**

Our work spans scales: from the regulatory mechanisms that establish and maintain
DNA methylation in single cells, to whole-genome methylomes across longitudinal
human cohorts of more than a thousand participants. We develop the computational
and statistical methods needed to interpret those data, and increasingly the
agentic systems needed to analyse them at scale.

**Affiliation line:**

Head of Epigenetics, Black Ochre Data Labs — The Kids Research Institute Australia
and the Australian National University. NHMRC Emerging Leadership Fellow. Adjunct
Senior Lecturer, The University of Western Australia.

> **[open]** Your edit changed ANU to UWA in the first sentence, which left UWA
> appearing twice, but you also gave me the ANU link. I've read that as a slip and
> restored ANU (linked to the JCSMR page you supplied), with UWA kept as the
> adjunct role. Correct me if you actually meant to drop ANU.

---

## 3. Hero figure — now a four-panel composite

**A · epigenome-wide scan.** Manhattan across the genome, one locus above the
significance threshold, labelled `chr8`. Zoom guides fan down to panel B.

**B · methylation at that locus.** Two groups, identical across flanking sequence,
separating at the promoter CpG island. Difference shaded and measured: `Δ 34%`.
Legend: `unaffected` / `affected`.

**C · variant association.** Local SNPs tested against methylation, plotted as
−log₁₀P on the same x scale as B. Lead mQTL marked with a diamond under the island.

**Annotation strip.** CpG density ticks, `CpG island`, transcript model, `8 kb`.

**D · variance partition.** Stacked bar: Genetics 34% · Cell composition 21% ·
Environment 18% · Development 12% · Unexplained 15%.

**Caption:**

A model of one part of the work, not real data. **A** An epigenome-wide scan
identifies a differentially methylated locus. **B** At that locus, two groups are
indistinguishable across the flanking sequence and separate at the promoter CpG
island. **C** Local genetic variants are tested for association with methylation;
a lead mQTL sits under the island. **D** The variance is partitioned between
genetics, cell composition, environment and development. Finding the difference is
the easy part; attributing it is the work.

> **[open]** The variance percentages are invented. Give me numbers that are
> plausible for your domain, or tell me to drop the numbers and label the segments
> only.
>
> **[open]** `chr8` is arbitrary — I moved it off chr1 because a peak at the far
> left edge made the zoom fan look wrong. Name a chromosome if you'd prefer.

---

## 4. Research questions

**Section eyebrow:** Open · **Heading:** Research questions
**Note:** What the current programme is pointed at.

---

**Q1 — How do genetic variants shape DNA methylation in health and disease?**

Most heritability for common disease lies outside coding sequence, in regulatory
variation we can detect but cannot yet interpret. Mapping the effects of variants
on the methylome is one route from association to mechanism — and the same logic
scales in both directions, from modelling a single patient's disease-causing
variant in cells to quantifying regulatory variation across a whole cohort.

**Q2 — How are environment and life history reflected in the epigenome?**

Exposure, development and social circumstance leave molecular signatures that
genotype alone cannot account for. Resolving those contributions, and establishing
which are stable, reversible or transmitted between generations, remains the
harder half of the problem.

> "marks" → "molecular signatures"; em dash removed as you asked. Title reworded
> to "How are environment and life history reflected in…" — your version had
> "How does the environment and life history become…" which disagrees in number.

**Q3 — Can epigenetic signatures become predictive and diagnostic tools that change clinical practice?**

An association that replicates across cohorts is not yet an assay. The distance
between the two is measurement precision, calibration in the population the test
will be used in, and evidence a clinical service can act on. We are interested in
signatures that survive all three.

**Q4 — How can agentic engineering pioneer the next revolution in genomics and precision medicine?**

Genomic analysis has always been rate-limited by the analyst, not the sequencer.
Agentic systems change that constraint — but only if what they produce is
verifiable. We are building and testing workflows where an agent specifies an
analysis, evaluates it against held-out data and negative controls, and emits a
record another researcher can reproduce and contest.

**Q5 — Can we build sovereign data systems for population-scale genomics?**

Analysis at this scale conventionally means moving data to compute. The
alternative is to move the analysis instead: federated systems in which data
remain under the control of the communities and institutions that hold them, and
participant privacy is preserved by the design of the system rather than resting
on it alone.

**Q6 — Can ethics and governance be enforced by the systems that hold the data, rather than relying on trust?**

Consent, access conditions and withdrawal are established through agreements, and
those agreements matter. But they are enforced today largely by good faith and
audit after the fact. We are interested in infrastructure where those conditions
are machine-readable and enforced at the point of access — a question that now
sits squarely at the intersection of genomics, privacy, data security and AI.

> Q5 and Q6 kept separate rather than merged: Q5 is the architecture, Q6 is what
> the architecture enforces. Q6's opening clause is deliberately non-dismissive of
> agreements, per your note. Say if you'd still rather merge them into one.

---

## 5. Current projects

**Section eyebrow:** In progress · **Heading:** Current projects
**Note:** Active, not complete.

---

**PROPHECY methylomes**

Whole-genome sequencing and EM-seq methylome profiling across a longitudinal
cohort of more than 1,200 participants with deep phenotyping, sampled repeatedly
over follow-up.

*Meta:* NHMRC Investigator Grant · 2025–2029

---

**Biomarkers for clinical use**

Developing epigenetic biomarkers for cardiometabolic risk, and extending polygenic
risk scores with methylation and quantitative trait (mQTL) information to improve
prediction in the populations where the burden is greatest.

*Meta:* MRFF Genomics Health Futures Mission · 2022–2024

> **[open]** You said this grant is attached to PROPHECY and to biomarker
> development, and asked for years. I've used 2022–2024 from your CV. If it has
> been extended, or if it should be listed against both items, tell me.

---

**Sovereignty by design**

Computational systems in which governance conditions agreed with participating
communities are enforced by the architecture itself, so that population-scale
analysis requires as little trust as possible from anyone involved.

*Meta:* Ongoing

---

**Verifiable agentic workflows**

Engineering the next generation of bioinformatics workflows for genomics and
precision medicine: agent-driven analyses built to be robust and verifiable by
construction, at a frontier that is moving faster than the methods literature can
document it.

*Meta:* Open, ongoing

> On the repetition: I cut it from the Method section rather than here, since you
> wanted both Q4 and this item strengthened. Method now leads on provenance and
> forking paths, and mentions agents only in passing.

---

## 6. Method

**Section eyebrow:** Method · **Heading:** Analysis you can audit
**Note:** How the computational work is structured.

---

A population-scale methylome study is a long sequence of analytical decisions:
coverage and quality thresholds, correction for batch and cell-type composition,
the choice of covariates, the handling of relatedness and population structure.
Each is defensible in isolation; together they define a garden of forking paths,
and a study can arrive at a confident wrong answer without any step visibly
failing.

We therefore treat provenance as a first-order requirement rather than
documentation added afterwards. Analyses record the decisions they made and the
evidence for them; assumptions are tested against negative controls and
permutation rather than asserted; and results carry enough of their own history
for another group to reproduce or contest them. Code is released openly wherever
the governing data agreements permit.

**Pull quote:**

The question is not whether a model can write the analysis. It is whether the
analysis can be shown to be wrong.

**Closing paragraph:**

The same requirement follows from the governance model. Where analysis must run
inside a sovereign environment and the underlying data cannot leave it,
auditability is not a matter of good practice — it is the only basis on which
anyone outside that environment can evaluate a result.

**Diagram:**

Cohort → Sequencing → Methylome → Model → Translation, with an agent layer beneath
containing three stages: `specify → evaluate → register`, annotated
"held-out data · negative controls · permutation · provenance for every call".
Methylome feeds into the agent layer; the agent layer returns into Model.

> **[open]** The three agent stages are my guess at your actual loop. Rename them
> to whatever you really do.

---

## 7. Selected publications

Eyebrow and section note both removed, as you asked. Heading is now just
"Selected publications".

---

Human iPS cells retain epigenetic memory of their somatic tissue of origin.
Transient passage through a naive state erases that memory and restores
developmental potential, correcting the cells both functionally and
epigenetically.

*Buckberry, Liu, Poppe, Tan et al. · Nature 620, 863–872 · 2023 · Paper · Patent*

---

A reprogramming roadmap resolving the transcriptional and epigenomic trajectories
of human somatic cell reprogramming, and identifying the route to induced
trophoblast stem cells.

*Liu et al. · Nature 586, 101–107 · 2020 · Paper*

---

Chromatin accessibility and transcription factor occupancy are reconfigured in
both transient and permanent modes during reprogramming, distinguishing the
changes that drive cell-fate conversion from those that accompany it.

*Knaupp & Buckberry et al. · Cell Stem Cell 21, 1–12 · 2017 · Paper · co-first author*

---

Gene regulatory dynamics of the human prefrontal cortex from gestation to
adulthood, resolved at single-cell resolution across the longest developmental
window in the brain.

*Herring et al. · Cell 185, 4428–4447 · 2022 · Paper*

---

Targeted methylation of hundreds of promoters shows that the transcriptional
response to promoter DNA methylation, and the stability of the mark itself, are
strongly context dependent.

*de Mendoza et al. · Genome Biology 23, 163 · 2022 · Paper*

---

Where epigenetic biomarkers for type 2 diabetes currently stand, and what is
required for them to be useful in global and Indigenous health rather than only in
the cohorts they were derived from.

*Munns, Brown & Buckberry · Frontiers in Molecular Biosciences 12, 1502640 · 2025 · Paper*

---

Consent and access management for genomic data must be revocable, machine-readable
and auditable to remain meaningful at scale; a review of how far current
technology supports that.

*Oliva et al. · GigaScience 13, giae021 · 2024 · Paper*

---

**Footer link row:** All publications · CV

> The TNT patent is now a second link on the Nature 2023 entry, pointing at the
> WIPO record (WO2021102500).
>
> **[open]** That's seven entries. Happy to cut back to five if it's running long
> — my candidates to drop would be Oliva and Herring.

---

## 8. Contact

**Label:** Get in touch · **Email:** sam.buckberry@thekids.org.au

**Note:**

Enquiries welcome from prospective PhD students and postdocs, and from anyone
interested in collaborating.

**Label:** Elsewhere · **Links:** GitHub · Scholar · ORCID · LinkedIn

**Portrait:** small, greyscale, in this block.

---

## 9. Colophon

Removed, as you asked.

---

## Your question about GitHub

Short answer: **don't mass-delete, do curate the front page.**

Deleting old repos is the one irreversible move here, and it costs you things you
can't get back — stars, issue history, inbound links from papers and Stack
Overflow answers, and the archaeology that shows you've been writing code for a
decade. A stale repo is not embarrassing. An *unexplained* stale repo is.

What I'd actually do, in order:

1. **Pin six.** The pinned row is the only part most visitors read. Pin your best
   maintained work, ideally including something tied to a paper. This is 90% of
   the benefit and takes ten minutes.
2. **Archive rather than delete.** GitHub's archive flag marks a repo read-only
   and greys it out. It reads as "finished", not "abandoned". Archive anything you
   won't touch again.
3. **Add a one-line description and a README to anything visible.** An unloved repo
   with a clear "code for Buckberry et al. 2023, unmaintained" README reads as
   scholarship. The same repo with no description reads as clutter.
4. **Hide the forks.** Stale forks are the main thing making a profile look messy,
   and the profile has a filter for exactly this. Delete forks you never committed
   to — those genuinely carry nothing.
5. **Write a profile README.** A short one naming what you work on and what's worth
   looking at. This is what makes a profile look deliberate.

On the private repos: releasing more is the right instinct, and for the
positioning you described it matters more than tidying. One well-documented public
repo attached to a paper does more than fifty archived ones. But don't let the
cleanup block the release — do the pinning and the profile README now, and release
as things become releasable.

The one thing I'd check before linking: whether anything public contains data,
paths or identifiers it shouldn't. If you want, I can audit the public repos for
that specifically.

---

## Applied from your last pass

- H1 unchanged; opening statement rewritten with development and cell identity
- Lead paragraph rewritten to span single cells to cohorts; "registers/declares" gone
- ANU link added; UWA spelled out
- Figure rebuilt as an A–D composite with mQTL and Manhattan panels
- Questions reordered — agentic engineering up to Q4, ethics to Q6
- Q2 title fixed, "marks" replaced, em dash removed
- Q3 rewritten out of comms register
- Q5 reflects federated/privacy-preserving principles without naming GA4GH
- "What's on the bench" → "Current projects"
- PROPHECY: >1,200, WGS, EM-seq, deep phenotyping; reference claim removed
- Biomarkers: PRS extended with methylation and mQTL; years added
- "by construction" → "by design"; agreements no longer dismissed; collaborator line removed
- Agent item levelled up to verifiable workflows
- Method rewritten in academic register; agent repetition cut here
- Diagram rebuilt with real stages instead of a single label
- Publications: Liu, Knaupp, de Mendoza (Genome Biology), Munns added; Khurana and
  the sponge paper removed; patent linked
- Contact note shortened; colophon dropped

---

## Not yet drafted

CV and publications pages are still untouched. Proposed change is restyling only —
funding amounts and journal names lose the accent colour so neither page reads as
a ledger. Say if you want either rewritten instead.
