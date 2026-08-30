# Homepage copy — round 6

Your edits are applied and live at the prototype. This file is the current state
of the copy. Same convention: each heading is a **slot** on the page, the text
under it is what renders, and **[open]** marks something I still need from you.

Resolved from your last pass is listed at the bottom so you can check I didn't
miss anything.

---

## 1. Header

**Site mark (top left, small):** Sam Buckberry

**Nav links:** Questions · Projects · Method · CV · Publications

> CV and Publications are now real page links (`/cv/`, `/publications/`), since
> dropping the publications section removed the only route to those pages. They
> also appear in the contact block at the foot, because the header scrolls away.
>
> **[open]** "Sam Buckberry" now appears twice near the top — once as the site mark,
> once opening the hero. Common enough on personal sites, and the two are set in
> different faces, but say the word and I'll drop the mark on the homepage only
> (keeping it on CV and publications, where it is the way back).

---

## 2. Hero

**Identity block — now the first thing in the hero, above the heading:**

Sam Buckberry, BSc, BHlthSc (Hons), PhD 

Head of Epigenetics, Black Ochre Data Labs — The Kids Research Institute Australia.
NHMRC Emerging Leadership Fellow. Adjunct Senior Lecturer at The University of
Western Australia and the Australian National University.

*Set small: the name in semibold sans at 0.95rem, the affiliation muted beneath it.
The heading below is still the largest thing on the page.*

**H1:**

Epigenetics of Health and Disease

**Statement (large serif, under the H1):**

The epigenome records the combined action of genetic variation, development, cell
identity and environment. We work to separate those contributions, and to use what
they reveal to understand, predict and prevent disease.

**Lead paragraph:**

Our work spans scales: from the regulatory mechanisms that establish and maintain
DNA methylation in single cells, to whole-genome methylomes across longitudinal
human cohorts. We develop the computational and statistical methods needed to
interpret those data, and increasingly the agentic systems needed to analyse them
at scale.

---

## 3. Hero figure — four-panel composite

**A · epigenome-wide scan.** Manhattan, one locus above threshold, labelled `chr8`.
Zoom guides fan down to panel B.

**B · methylation at that locus.** Two groups, identical across flanking sequence,
separating where CG density is highest. `Δ 34%` retained — the one hard number left.

**C · variant association.** −log₁₀P on the same x scale. Numeric axis removed.
Lead mQTL sits ~6 kb distal, outside the shaded region.

**CG density track.** Labelled at top left, like the other panels. `CpG island` text
removed — the band plus the tick density carries it without saying it. A sashimi
arc runs from the variant to the target CpG, annotated only `6 kb`.

**D · variance partition.** Proportional segments, no percentages: Genetics ·
Cell composition · Environment · Development · Unexplained.

**Caption:**

Model data of one part of the work. Epigenome-wide scans
identify differentially methylated regions, and testing nearby genetic variants for association with methylation puts the
lead mQTL several kilobases away, acting on the CpG at the centre of the change
from a distance. This begins the process of partitioning variance between genetics, cell-type
composition, environment and development.
## 4. Research questions

**Section eyebrow:** Open · **Heading:** Research questions

---

**Q1 — How do genetic variants shape DNA methylation in health and disease?**

Most heritability for common disease lies outside coding sequence, in regulatory
variation we can detect but cannot yet interpret. This is also the case for many rare and undiagnosed diseases. Mapping the effects of variants
on the methylome is one route from association to mechanism — and the same logic
scales in both directions, from modelling a single patient's disease-causing
variant in cells to quantifying regulatory variation across a whole cohort.

**Q2 — How are environment and life history reflected in the epigenome?**

Exposure, development and environment leave molecular signatures that
genotype alone cannot account for. Resolving those contributions is the harder half
of the problem, and it turns on a distinction that is easy to state and difficult
to establish: which signatures are stable, which are reversible, and which merely
record an exposure rather than mediating its effect on risk.

> Transgenerational inheritance is gone. The replacement ends on the
> mediator-versus-marker problem, which is a live methodological question in
> exposure epigenetics and carries none of the same baggage. I checked the whole
> page: no remaining mention of inheritance across generations.

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

**Q5 — Can sovereign data systems enforce the governance they run under?**

Analysis at population scale conventionally means moving data to compute. The
alternative is to move the analysis instead, so that data remain under the control
of the communities and institutions holding them. That only helps if the governance
travels with it. Consent, access conditions and withdrawal are set out in agreements
that matter, but are enforced today largely by good faith and audit after the fact.
Making them machine-readable and enforced at the point of access is a problem
sitting squarely at the intersection of genomics, privacy, data security and AI.

> The two governance questions are merged, so the list is five rather than six. The
> merge keeps the federated-architecture opening and the enforcement problem, and
> drops the seam between them. The Method section still carries the distinct point
> about bounding an analysis to its approved question.

---

## 5. Current projects

**Section eyebrow:** In progress · **Heading:** Current projects
**Note:** Active, not complete.

---

**PROPHECY epigenetics program**

PROPHECY is an Aboriginal longitudinal cohort established to investigate
cardiometabolic disease: multi-omic, deeply phenotyped, and now more than 1,200
participants. Our program within it is whole-genome sequencing and EM-seq methylome
profiling, sampled repeatedly across follow-up.

> Placeholder framing — you said you'd add the detail. Correct anything wrong about
> the cohort description before it goes live.

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

**Data sovereignty by design**

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

**Section eyebrow:** Method · **Heading:** Reproducibility without the data
**Note:** Why the usual guarantee is unavailable, and what replaces it.

---

The standard guarantee of reproducible research is that another group can obtain
the data and re-run the analysis. For human genomic data that guarantee rarely
holds. Access is restricted by construction — by consent conditions, ethics
approvals and data-sharing agreements — and as genetic privacy becomes a sharper
concern, those restrictions are tightening rather than easing.

Where analysis must run inside a controlled environment and the underlying data
cannot leave it, auditability is not a matter of good practice — it is the only
basis on which anyone outside that environment can evaluate a result.

Meeting that standard needs tooling the field does not yet have. The documents that
govern this work — approvals, consent conditions, governance agreements — are
written in prose and read by people. The analyses they govern run in code and are
read by almost no one. Closing that gap means translating governance into
constraints a system can enforce, so that an analysis is bounded by the question it
was approved for and cannot quietly extend into questions it was not.

Within those bounds, the burden falls back on the analysis itself. A
population-scale methylome study is a long sequence of decisions: coverage and
quality thresholds, correction for batch and cell-type composition, the choice of
covariates, the handling of relatedness. Each is defensible alone, and together
they define a garden of forking paths in which a study can arrive at a confident
wrong answer without any step visibly failing. What can leave the environment is
the record — which path was taken, what was tested against negative controls and
permutation, and what those tests ruled out.

**Pull quote:**

The question is not whether a model can write the analysis. It is whether the
analysis can be shown to be wrong.

**Diagram:** removed.

> The framing is now field-level rather than local: the problem belongs to human
> genomic data generally, not to your cohorts specifically, and genetic privacy
> tightening is named as the direction of travel.
>
> The third paragraph is the new claim, and it is the one that distinguishes this
> section from Q6. Q6 asks who can get *in*; this asks what you are permitted to
> *ask* once you are there. Stated as a gap in the field's tooling rather than as
> something already solved.

## 7. Selected publications — removed

Section dropped for a cleaner, more minimal landing page. The publications page is
untouched and still lives at `/publications/`, linked from the nav and the contact
block.

---

## 8. Contact

**Label:** Get in touch · **Email:** sam.buckberry@thekids.org.au

**Note:**

Enquiries welcome from prospective PhD students and postdocs, and from anyone
interested in collaborating.

**Label:** Elsewhere · **Links:** CV · Publications · GitHub · Scholar · ORCID · LinkedIn

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

## Applied from the previous pass

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
