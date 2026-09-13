# get-RPF read-structure inference — specification (v1)

- **Status:** draft for implementation, 2026-09-13
- **Supersedes:** `docs/rpf_extraction_rewire_implementation_guide.md`, `docs/structure_inference_validation.md` (both preserved on branch `wip/codex-read-structure-rework`)
- **Background:** `docs/riboseq_library_observability_reference.md` (same branch); §12 lists corrections to it

## 1. Goal

Given a Ribo-seq FASTQ, decide which bases of each read are the biological fragment (the RPF) and which are technical, and be able to say why. Emit a per-read transform only when the evidence supports one. Otherwise, report exactly which question is unresolved.

The product is an **architecture**: an ordered list of blocks with lengths, each carrying its evidence. Trimming follows from the architecture; it is never a separate heuristic.

**Principle:** every analysis answers a named question. A more expensive analysis runs only if a cheaper one left that question unresolved *and* the answer would change the emitted bases.

### Not in v1

- Paired-end and index-read inference: detect and report only.
- Per-read assignment in mixed libraries: detect and report only.
- Naming the wet-lab protocol. Templates are hypotheses to test, not answers.
- Fragment policies beyond the named default (`monosome_20_40`).

## 2. The read model

A read is a straight walk through the library molecule. The order of elements is fixed by the chemistry; only their lengths vary.

```text
R1:  [5' technical] [ insert (RPF) ] [3' technical] [adapter] [index] [P7] [no signal]
     |- read-start frame -|          |------- adapter frame ------|
```

The implementation relies on four rules:

1. **Two coordinate frames.** 5′ technical blocks sit at fixed cycles from the read start. 3′ technical blocks sit at fixed offsets *before the adapter start*. In raw reads the read end has no structural meaning, because it lands in the adapter, index, P7 or poly-G. Analyse 3′ structure relative to the read end only for reads that were already trimmed.
2. **Variable bases occur at junctions.** UMIs, randomised ligation bases and non-templated additions (NTAs) sit only in two places: between the read start and the insert, or between the insert and the constant part of the adapter. Search nowhere else.
3. **Arithmetic.** `adapter_start = len(5' technical) + len(insert) + len(3' technical)`. Subtracting a typical footprint length from the most common adapter start estimates how much technical sequence to look for. This only sizes the search window. It never places a boundary, so expected RPF length never defines the RPF.
4. **Circularisation splits randomisation.** In circularisation protocols the RT primer's random bases appear at the read start and the linker's random bases appear 3′ of the insert. A "distributed UMI" is simply two ordinary blocks, one in each frame.

Block types:

| Type | Frame | Across the library | Within copies of one fragment\* |
|---|---|---|---|
| `umi` / `random` | either | even base mix, high diversity | low (≈0.25) |
| `fixed` (spacer, linker constant, template-switch motif) | either | one sequence | high |
| `barcode` (inline) | either | a few discrete values | high |
| `nta` | at a junction | skewed base mix, 0–3 nt, ragged length | low |
| `tail` (poly(A) etc.) | adapter | homopolymer, variable length | — (identified by position) |
| `adapter` | defines the adapter frame | constant; absolute position varies | high |
| `insert` | between the frames | diverse | high |

\*Measured by the pileup (§5.3c). Only the two columns together separate biology from technical sequence. `insert` and `fixed` both agree within copies of a fragment, but only `insert` differs between fragments.

## 3. The questions

Inference is organised around five questions. Each gets its own answer, status, evidence and rejected alternatives (§4.1), and the architecture is their combination.

| ID | Question | Changes emitted bases? |
|---|---|---|
| Q1 | Is there a 3′ anchor, what is it, and where does it start in each read? | yes |
| Q2 | What sits between the read start and the insert (type, length)? | yes |
| Q3 | What sits between the insert and the anchor (type, length)? | yes |
| Q4 | Where are the UMIs, if any, so they can be kept? | no (metadata) |
| Q5 | Does the emitted insert behave like Ribo-seq, and which fragment class is it? | no (validation) |

- **Escalation rule:** run a more expensive method only while Q1, Q2 or Q3 is unresolved.
- **Stop rule:** stop as soon as further work cannot change the transform.

## 4. Statuses

Every answer carries one status. Where the answer is numeric, it also carries a confidence (0–1) and the underlying measurement.

| Status | Meaning |
|---|---|
| `resolved` | One answer is supported, with an exact value |
| `interval` | Supported, but only to a range (e.g. insert end 29–31 at a poly(A) junction) |
| `ambiguous` | Evidence supports two or more incompatible answers |
| `conflicting` | Methods disagree (e.g. pileup vs alignment) |
| `underpowered` | The method applies but the sample had too little data (e.g. too few pileup groups) |
| `not_observable` | The information is not in the supplied reads (removed before deposit, or the UMI is in I1) |

`underpowered` and `not_observable` never turn into "absent". A UMI is reported `absent` only when all of the following hold:
- Q2 and Q3 are both `resolved` with no random block.
- Step 0 found no header or index-read UMI channel.

### 4.1 Provenance

A status alone is not an answer. Every answer records the evidence that produced it and the alternatives it considered, each with the reason it was rejected. This record, not a confidence number, is what makes a boundary defensible.

```text
Q2  5' technical block = 4 nt (random)                 status: resolved
    evidence
      - positions 1-4: even base mix across the library (entropy 1.98 bits)
      - positions 1-4: agreement within repeated-fragment groups 0.26-0.28
      - agreement rises to 0.96 at position 5 (612 groups)
      - no fixed motif at positions 1-6
    alternatives
      - 0 nt: contradicted - positions 1-4 disagree within groups
      - 5 nt: contradicted - position 5 agreement 0.96 is at biology level
      - 3 nt: contradicted - position 4 agreement 0.27 is at random level
```

## 5. Pipeline

Steps 0–5 run on a bounded, deterministic sample (default: the first 300k reads, configurable). The architecture is shared by every read, so no inference step needs the whole library. Only the transform (§7) reads the whole file.

```text
0 Observe → 1 Profile → 2 Anchor (Q1) → 3 Blocks + pileup (Q2–Q4) → [4 Align, only if needed]
          → 5 Assemble → 6 Transform → 7 Validate (Q5)
```

### 5.0 Observe (free)

Record which information channels exist:
- R1, R2, I1, I2
- a UMI field in the read name
- cycle count
- instrument chemistry (two-colour instruments give terminal poly-G)
- any supplied metadata or seqspec

Metadata is stored as a *claim* next to the observed evidence, never as truth. Anything with no channel gets `not_observable`.

**Output:** `Observation`.

### 5.1 Profile (one pass)

Measure on the sample:
- length distribution and fraction of reads at the modal length
- per-cycle base mix and entropy, and per-cycle quality
- duplicate fraction
- top k-mers and how much their positions vary
- homopolymer runs

**Output:** input state `raw_fixed_length` | `trimmed` | `mixed`, with the evidence. Profiling never trims.

### 5.2 Anchor — Q1 (one pass)

1. **Known adapters, found anywhere in the read.** Seed with k-mers, then extend while allowing mismatches (`max_mismatch_per_10nt`, default 1).
   - A **full** match may occur anywhere at or after `min_insert` (default 15).
   - A **partial** match is accepted only where it runs to the read end, with at least `min_partial_overlap` bases (default 10).
   - A match starting before `min_insert` is counted as a dimer, not an anchor.
   - Record per read: adapter id, start, matched length, mismatches.
2. **Chain check.** Where reads extend past the adapter, confirm the downstream order: adapter → index → P7 → poly-G. Catalogue entries that are downstream elements are marked `role: downstream` and can never be the anchor. `ATCTCGTATGCC…` (P7) is one of these.
3. **De novo discovery.** Use this if no known adapter reaches `min_anchor_support` (default 0.3). Look for an over-represented k-mer whose absolute position varies across reads, whose occurrences are always followed by the same continuation, and whose upstream sequence is diverse.
4. **Several strong anchors.** Record each with its read fraction and set Q1 to `ambiguous` (possible mixed library). v1 stops at reporting.
5. **Tails.** A homopolymer run directly upstream of the adapter is a `tail`. Its start is a second, per-read anchor, and the insert's 3′ end is measured from it.

**Output:** `Anchor`: sequence, source, per-read start distribution, support fraction.
- Trimmed input with no anchor → Q1 `not_observable`. This is not a failure.
- Raw input with no anchor → Q1 `underpowered` or `ambiguous`, as the evidence dictates.

### 5.3 Blocks and pileup — Q2, Q3, Q4 (seconds)

**(a) Search windows.** Estimate technical length as `T ≈ mode(anchor start) − typical_footprint` (config, default 30; per organism if known). Search from 0 to `max(T + 4, 16)` nt, both at the read start and before the anchor. For trimmed input, search the same widths at both read ends.

**(b) Profile both frames.** For each position in each window, measure the base mix, entropy, and number of distinct 4–6-mers. Classify runs of positions as `fixed`, `few_values`, `random`, `homopolymer`, `skewed` or `diverse`.

**(c) Pileup (the library as its own reference).** Always run this on the bounded inference sample, never on the whole library. At that size it is cheaper to run than to decide whether to run: the prototype took 2–3 s for 300k reads in pure Python and placed 40% of SRR1944950 reads into about 560 groups.

1. Count k-mers (default k = 14) in the region upstream of the anchor. Adapter k-mers must never become seeds. Drop low-complexity seeds.
2. Assign each read to its highest-ranked seed. Keep groups with at least `min_group_distinct` distinct sequences (default 5). Collapse identical sequences first, because PCR duplicates carry no information.
3. Line up each group's members on the seed. For every read position, in both the read-start frame and the adapter frame, score whether the base equals the leave-one-out consensus of the other distinct sequences covering that offset. Require at least 3 others.
4. Report agreement per position, the base mix of the disagreeing bases, and the number of groups.
5. If there are fewer than `min_groups` informative groups (default 100), the pileup is `underpowered`.

**(d) Classify and place boundaries.** Combine (b) and (c) using the §2 table.
- A boundary is the change point where agreement moves between the biology level (`agree_bio`, default ≥ 0.9) and the random level (`agree_random`, default ≤ 0.4).
- A sharp change point is `resolved`.
- A change point spread over several positions is `interval`. NTA raggedness and A-ending inserts next to a tail both do this.
- Every candidate boundary position in the window is recorded as an alternative, with its verdict (§4.1).

**Reference results the implementation must reproduce** (SRR1944950 with spiked controls, prototype `pileup_probe.py`):

| Case | Agreement from read start | Agreement from read end |
|---|---|---|
| real reads | pos 1: **0.64**, then 0.93–0.98 | 0.94–0.98 |
| + 5 nt random 5′ | pos 1–5: 0.26–0.27, pos 6: 0.80, then 0.97+ | 0.97–0.99 |
| + 4 nt random 3′ | pos 1: 0.78, then 0.97+ | last 4: 0.26, then 0.96+ |
| + untemplated G in 70% of reads | pos 1: 0.45, pos 2: 0.74, then 0.93+ | 0.94+ |

In the real reads, the disagreeing first bases are 60% T. That makes a 5′ NTA the leading explanation in a large fraction of reads. It is a validation target (§10.2), not yet a trimming rule.

### 5.4 Align (escalation only)

**When to run:** only if both conditions hold:
- Q2 or Q3 is `underpowered`, `ambiguous` or `conflicting`;
- a reference was supplied.

**How:** one local alignment of at most 10k reads (STAR `--alignEndsType Local`, reusing `processors/alignment.STARAligner`). Read three things:
- the 5′ soft-clip length histogram
- the 3′ soft-clip length histogram
- first- and last-base mismatch rates

**Interpretation:**
- A discrete soft-clip mode at k supports a k-nt technical block.
- A spike of 1-nt clips or mismatches with a skewed base mix supports an NTA.
- Do not compare mapping rates across many candidate trims.
- If alignment and pileup disagree, set `conflicting` and keep both answers. Alignment never silently overrides the pileup.

### 5.5 Assemble

Build the architecture from the per-question answers.

Known templates enter as hypotheses. Each template block is checked against the observed blocks, and the template is marked:
- `supported`: used as is;
- `partially_supported`: completed from the observations;
- `rejected`.

Where Q2 or Q3 has two surviving answers, emit both architectures with status `ambiguous`. The leading one may be displayed but is never applied.

There is no global additive score. Transform confidence is bounded by the least-resolved question that affects emitted bases (Q1–Q3). A well-supported adapter cannot make up for an unresolved 5′ boundary.

### 5.6 Transform

See §7. This is the only step that removes bases.

### 5.7 Validate — Q5

On the emitted inserts, measure:
- the length distribution against the named fragment policy
- the dimer fraction (adapter or tail directly after the 5′ block)
- the low-complexity fraction, reusing `release.py` and `identity_screen`
- with alignment: rRNA fraction and frame preference by length

**Output:** an identity call: `riboseq_likely`, `riboseq_possible`, `rna_like`, `technical_failure` or `insufficient_evidence`.

Validation never moves a boundary. A validation failure is reported against the architecture that produced it.

## 6. Poly(A)

Poly(A) is a structural clue. Classify it by position:

| Pattern | Interpretation |
|---|---|
| short insert → A-run → adapter | tailing library (D-Plex, SMARTer smRNA, CATS); the A-run is a `tail` block |
| 5′ technical → A-run → adapter, no insert | tail dimer: drop those reads and count them |
| long inserts of every length, containing A-runs | RNA-like material: Q5 `rna_like` |
| terminal poly-G after the construct | empty cycles on a two-colour instrument, not a tail |

**Junction convention:** cut at the first base of the A-run. Inserts that genuinely end in A lose those bases. Report the boundary as `interval`, naming the convention.

## 7. Transform and outputs

### 7.1 Contract

The `Architecture` object is the contract between inference and extraction. A single function serves both audit and production:

```python
apply(architecture, read) -> Accepted(insert, qual, umi) | Rejected(reason)
```

- It runs per read, streaming.
- Audit computes the same results but does not write the RPF FASTQ.
- Both audit and production record a transform hash over the architecture plus the fragment policy.

### 7.2 When a transform is emitted

A transform is emitted only when Q1, Q2 and Q3 are each `resolved`, or `interval` under a named convention (§6).

Otherwise the architecture is reported and the transform is withheld. This uses an exit status distinct from an error.

NTAs are reported but not trimmed in v1, because whether a given read carries one can't be decided without a reference. The report states the estimated rate.

### 7.3 Outputs

| Output | Content |
|---|---|
| `<sample>.rpf.fastq.gz` | Accepted inserts within the fragment policy, with original read names and the qualities of the retained bases only |
| read-name UMI | Appended as `_<UMI>` (the umi_tools convention); 5′ then 3′ UMI, concatenated in read order |
| `<sample>.structure.json` | Observation, per-question answers with evidence and rejected alternatives, architecture(s), and rejection counts by reason |
| `<sample>.seqspec.yaml` | The architecture as seqspec: region types `umi`, `barcode`, `linker`, `cdna`, `poly_A`, `adapter`, with variable min/max lengths |
| `<sample>.collapsed.fa` | Optional, as today |

### 7.4 Explanations

Every emitted boundary carries a one-sentence reason generated from its evidence. For example:

> Insert starts at position 6: positions 1–5 are a fixed length from the read start, evenly mixed across the library, and disagree between copies of the same fragment (agreement 0.27 → 0.97 at position 6; 612 groups). Alignment was not run because the architecture was resolved without a reference.

## 8. Data model (sketch)

```python
class Status(str, Enum):
    RESOLVED = "resolved"; INTERVAL = "interval"; AMBIGUOUS = "ambiguous"
    CONFLICTING = "conflicting"; UNDERPOWERED = "underpowered"; NOT_OBSERVABLE = "not_observable"

Frame = Literal["read_start", "anchor"]   # anchor = adapter start, or the read end for trimmed reads

@dataclass(frozen=True)
class Evidence:
    source: str      # profile | anchor | blocks | pileup | alignment | template | metadata
    metric: str      # e.g. "agreement_at_position"
    value: Any       # the continuous measurement, always reported
    n: int           # reads or groups behind it
    note: str = ""

@dataclass(frozen=True)
class Alternative:
    value: Any       # the rejected answer, e.g. a boundary position
    verdict: Literal["contradicted", "weakly_compatible", "untested"]
    reason: str

@dataclass(frozen=True)
class Answer:
    question: Literal["Q1", "Q2", "Q3", "Q4", "Q5"]
    value: Any
    status: Status
    evidence: tuple[Evidence, ...]
    alternatives: tuple[Alternative, ...]
    explanation: str                 # §7.4

@dataclass(frozen=True)
class Block:
    type: Literal["umi", "random", "fixed", "barcode", "nta", "tail", "insert", "adapter"]
    frame: Frame
    length: tuple[int, int]          # (min, max); tail and insert may vary
    sequence: str | None             # fixed and adapter blocks only
    keep_as_umi: bool
    status: Status
    evidence: tuple[Evidence, ...]

@dataclass(frozen=True)
class Architecture:
    blocks: tuple[Block, ...]        # in read order
    source: Literal["inferred", "template", "template_completed"]
    status: Status
    fragment_policy: str

@dataclass
class StructureReport:
    observation: Observation
    answers: dict[str, Answer]       # Q1..Q5
    architectures: list[Architecture]   # more than one only when ambiguous
    transform_emitted: bool
    counts: dict[str, int]           # input, accepted, rejected by reason
```

All thresholds live in one `InferenceConfig` dataclass. Each has a name, a default and a docstring, and the report prints the measured value next to every thresholded decision.

## 9. Code layout

Build a new package, and keep the existing `extract` path working until phase P6.

```text
src/getRPF/core/structure/
  model.py        Status, Evidence, Alternative, Answer, Block, Architecture, StructureReport
  config.py       InferenceConfig
  observe.py      §5.0–5.1
  anchors.py      §5.2
  blocks.py       §5.3 a, b, d
  pileup.py       §5.3 c
  align.py        §5.4 (wraps processors/alignment.STARAligner)
  assemble.py     §5.5, including template testing
  transform.py    §7.1 apply()
  seqspec_io.py   Architecture <-> seqspec YAML
  report.py       structure.json and explanations
tests/test_structure/
  sim.py          synthetic read generator with known truth
  test_*.py
```

**CLI:**
- `getRPF infer-structure <fastq> -o <dir>` writes the report and seqspec, but no FASTQ.
- `getRPF extract --architecture <seqspec>` applies a resolved architecture.
- Existing `extract` behaviour is unchanged until P6.

**Why not reuse `ReadArchitecture`:**
- It stores UMIs only as absolute read coordinates.
- `seqspec_loader` turns any UMI or barcode after the insert into `post_rpf_trim_bases`, so the UMI is discarded rather than kept.
- It has no tail, NTA or variable-length blocks.

Seqspec remains the serialisation and interchange format; it does not constrain the internal model. The template YAMLs in `src/getRPF/architectures/` stay, and `seqspec_io.py` reads them into `Architecture`.

**Not carried over** from the WIP branch or the current extractor:
- read-end-only adapter matching (`RPFExtractor._find_adapter_prefix`)
- architecture selection by adapter-hit fraction
- HMM segmentation as a source of boundaries
- the 17-alignment, mapping-rate trim search (`_bounded_alignment_trim_search`)
- additive `EvidenceSummary` scoring
- rewriting the FASTQ from collapsed counts, which loses read names and qualities

## 10. Testing

### 10.1 Synthetic tests (network-free, run in CI)

`sim.py` builds a random transcriptome, samples footprints with skewed (Zipf) abundance so fragments recur as they do in real Ribo-seq, wraps them in a declared architecture, and simulates sequencing at a chosen read length.

Tests assert on blocks, boundaries **and statuses**:

| # | Architecture | Expected |
|---|---|---|
| 1 | insert · adapter | no technical blocks; insert exact |
| 2 | UMI5 · insert · adapter | 5′ `umi` 5, resolved |
| 3 | insert · random4 · adapter | 3′ `random` 4, resolved |
| 4 | UMI5 · insert · random4 · adapter | both blocks |
| 5 | UMI12 · fixed4 · insert · A(5–30) · adapter, 100 nt reads | `tail` block; insert end `interval` under the convention |
| 6 | NTA(0–2, skewed) · insert · adapter | `nta` reported, not trimmed |
| 7 | insert · NTA(1) · adapter | 3′ `nta` |
| 8 | insert · random5 · barcode5 (one value) · adapter | 3′ `random` kept as UMI, then `fixed`/`barcode` |
| 9 | as 8, with 4 barcode values | `barcode` with `few_values` |
| 10 | already-trimmed inserts | Q1 `not_observable`; Q2/Q3 from the pileup |
| 11 | adapter not in the catalogue | de novo anchor |
| 12 | two architectures, 70/30 | Q1 `ambiguous`; transform withheld |
| 13 | very low duplication | pileup `underpowered`; UMI never called absent |
| 14 | insert ending in AAA, then a tail | `interval`, never exact |
| 15 | reads shorter than the insert | no 3′ cut invented |
| 16 | insert · adapter · index · P7 in 100 nt reads | anchor at the adapter, never at P7 |

A test fails when the result is **more certain than the truth allows** (for example, an exact boundary in #14 or UMI `absent` in #13), not only when it is wrong.

### 10.2 Real panel

**Gold set rules:**
- Scope each claim to the run, sample or study it actually came from, and cite the source.
- A blank cell means unknown, never `False`.
- Keep three separate truth fields: protocol claim, read-observable structure, expected RPF output.
- A value the project measured itself is marked as a measurement, with the date and method, and stays separate from any published claim.

**Corrections carried into the gold set:**
- **SRR1944950:** the WIP gold quoted McGlincy & Ingolia 2017, but the run belongs to GSE67387 (Nedialkova & Leidel 2015) and was adapter-trimmed before deposit. Its protocol UMI claim is `unknown`; the observed structure is "trimmed, candidate 5′ NTA".
- **SRR3945920:** this is `DualLigation_RibosomeProfiling_rep1`, not a 4+3N sample. The 4+3N case is SRR3945930 (`4+3N_RibosomeProfiling_rep1`).
- **SRR23242345:** ENA lists it as single-end (101 nt), so it belongs in the single-end panel.

**Panel v1:**

| Run | Why it's in the panel |
|---|---|
| SRR1944950 | trimmed input; candidate 5′ NTA |
| SRR3945920 | 50 nt reads; adapter runs off the read end |
| SRR3945930 | Lecanda 4+3N: circularisation, random blocks in both frames |
| SRR12693498 | D-Plex, 100 nt; 5′ UMI plus tail |
| SRR23242345 | OTTR, 101 nt; 5′ UMI, adapter mid-read |
| McGlincy-style run (to be sourced) | 3′ UMI plus inline barcode |

Confirm the SRR1944950 NTA by alignment (§5.4) before any trimming behaviour depends on it.

## 11. Delivery plan

Each phase lands on `dev` through a feature branch, with its tests, and passes `make lint`, `make typecheck` and `make test`.

**Order:** P0 → G → P1 (minimal) → **P3 gate** → P2 → P4 → P5 → P6 → P7 → P8.

P3 comes before P2 on purpose. Its tests take anchor positions from the simulation truth, so the gate tests the pileup alone. The pileup is the novel part; prove or falsify it before building the rest around it. No further framework design until the gate is decided.

| Phase | Scope | Done when |
|---|---|---|
| P0 | Preserve the uncommitted Codex rework on `wip/codex-read-structure-rework` (excluding `search_cache/`, `search_results/` and HPC outputs over 5 MB), keep a patch + tarball backup outside the repo, reset `dev` to `origin/dev` | WIP commit verified against the backup; `dev` at `24cf600` |
| G | Gold set v1 per §10.2 | every claim has a scope and a source; measurements marked as such |
| P1 | Minimal `model.py`, `config.py`, `sim.py` | simulator produces reads with per-read truth |
| P3 | `pileup.py` plus the boundary calls it needs — **go/no-go gate** | cases 1–4, 6–9 and 13 pass, statuses and alternatives included; SRR1944950 reproduces the §5.3 table |
| P2 | `observe.py`, `anchors.py` (including tails) | cases 1, 5, 10, 11, 12, 15, 16 |
| P4 | `blocks.py`, `assemble.py`, `report.py`, `infer-structure` CLI | all 16 cases (5 and 14 need P2 tails) give the correct architecture and status, with explanations |
| P5 | `transform.py`, `seqspec_io.py`, `extract --architecture` | streaming output keeps read names and qualities; audit and production hashes match; templates round-trip through seqspec |
| P6 | Real panel with the gold set; decide whether `extract` defaults to the new path | panel report per §10.2 |
| P7 | `align.py` escalation | SRR1944950 NTA confirmed or refuted |
| P8 | R2 overlap; per-read assignment in mixed libraries | — |

If P3 fails its gate, meaning the pileup can't robustly recover the synthetic blocks, stop and revisit the design before building anything else.

## 12. Corrections to existing docs

`riboseq_library_observability_reference.md` (on the WIP branch):
- PMC5024760 is Lecanda et al. 2016, not "Kwak and colleagues".
- Its claim that a leading field is usually unobservable from the FASTQ alone is too strong. The pileup observes it whenever fragments recur.
- Its SRR1944950 paragraph rests on a protocol claim from the wrong study.

`rpf_extraction_rewire_implementation_guide.md` (on the WIP branch):
- Invariant 6 ("adapter trimming requires an observed end-boundary match") causes the long-read failures. §5.2 supersedes it.
