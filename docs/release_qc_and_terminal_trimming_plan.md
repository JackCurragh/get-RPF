# getRPF next major release: trustworthy pre-alignment RPF extraction

## Scope

Define the algorithm, checks, output contract, and milestones for a getRPF release whose job is to obtain clean ribosome-protected fragments (RPFs) from raw FASTQ with defensible, machine-readable evidence, and to refuse rather than guess when the evidence is insufficient. This stage is pre-alignment only. A second, alignment-based gate later in the workflow confirms biological identity; this document is explicit about which claims that split lets each stage actually make.

This supersedes the framing of the previous "release QC and terminal trimming plan." It keeps that plan's output-class and provenance goals but re-grounds the trimming logic in convergent evidence and an explicit refusal path, and it reconciles the milestones with capabilities already present in the codebase.

## 1. The epistemic contract

Every extraction implicitly asserts two independent claims:

1. Structural — the biological insert occupies read coordinates [i, j]; everything outside is adapter, UMI, barcode, or linker and can be removed.
2. Biological — that insert population is ribosome footprints, i.e. the library is Ribo-seq.

Pre-alignment FASTQ can earn the structural claim and only screen the biological one. The single signal that is specific to Ribo-seq — sub-codon triplet periodicity relative to start codons — is frame-relative and therefore alignment-only. Length, composition, contamination fraction, and complexity are all consistent-with signals shared by small-RNA, degradome, miRNA, and over-digested libraries.

The posture that follows is asymmetric and must be reflected in the code's outputs:

* Pre-alignment is strong at rejection: it can confidently say "the insert cannot be located" or "this is not a clean small-RNA/footprint library."
* Pre-alignment is weak at confirmation: the strongest positive verdict it may emit is structurally clean, footprint-length-consistent, contamination-consistent — biological identity pending the alignment gate.

The tool must never emit "confirmed Ribo-seq" from FASTQ alone. Every release-compatible sample carries `biological_confirmation: pending_alignment_gate`. Making the structural claim rigorous and the biological claim explicitly provisional is the design, not a shortcoming.

## 2. Design principles

Convergence, not any single estimator. The insert boundary is estimated by several independent methods (Section 6). Concordance among them is the trust signal; divergence is itself information — it means the architecture is not understood, and the sample is routed to review rather than averaged into a false decision.

5′ and 3′ risk are not symmetric. The 5′ boundary defines frame and P-site interpretation at the next gate; a 1-nt error there corrupts downstream periodicity. The 3′ boundary is merely where ligation occurred. The confidence bar for 5′ trims (UMI/barcode/linker removal) is therefore strictly higher than for 3′ adapter trims, and 5′ trimming is never performed on de novo terminal-base bias alone — it requires a named element from a seqspec or a high-posterior segment.

Confidence scales with support. A 0.998 terminal-base fraction over 280k reads is certain; the same fraction over a few hundred reads is noise. This binds especially hard on per-length rules, where rare length classes will otherwise manufacture spurious trims. Every rule carries an explicit support floor.

Refusal is a first-class output. `hold` and `needs_*` are not failures of the tool; they are correct verdicts. The performance target is to make `safe_to_trim` as broad as possible without ever trimming into the RPF body — recall is subordinate to that safety constraint.

Infer on a subsample, apply on a stream. All inference runs on a bounded subsample (~10⁵–10⁶ reads). Application streams over the full file using frozen rules. Inference and application are separate passes with separate cost profiles.

Keep the gates separate. No reference genome or transcriptome is consulted at this stage. Contamination screening uses reference k-mer sets (rRNA/tRNA/sno), not alignment. Anything requiring frame or coordinates is reserved for the alignment gate (Section 12).

## 3. Reconciliation with the current codebase

This plan builds on infrastructure that already exists and removes duplication rather than adding a third implementation.

Reuse:

* `SignalProcessor` (`signals.py`) already computes coverage-normalised per-position entropy and composition in both 5′- and 3′-anchored coordinates plus dinucleotides. This is the correct entropy engine and becomes the single source of per-position statistics.
* The HMM segmenter (`segmenter.py`, `plot-hmm`) is the basis for architecture segmentation; it needs a concordance check and a refusal path, not more emission features.
* `decide-trim` / `TrimDecider` (`consensus.py`) already reconcile architecture detection with alignment; the pattern generalises to multi-estimator concordance.
* `seqspec_generator.py` / `seqspec_loader.py` / `architectures/` provide protocol assignment and promotion.
* DuckDB ingestion provides the cohort-level machine-readable store.

Fix / retire:

* `check.py` normalises per-position frequencies by total reads, not by reads covering the position, which makes `InformationContentCheck` deflate entropy at variable-length tails and produce false `low_complexity` calls. Retire that bespoke entropy path; route cleanliness QC through `SignalProcessor`.
* `handle_extract_rpf`'s "Whole Shebang" mode computes a consensus trim and then discards it (`TODO: Pass override trims to extract_rpfs`). There is no per-length override-trim parameter threaded into extraction. This apply path is the real critical path for release provenance and is promoted to a first-class milestone.
* No per-length terminal diagnostics exist anywhere; `SignalProcessor` pools all lengths. This is the one genuinely new signal to build.

## 4. High-level algorithm

Stage 0 — Sketch. On a bounded subsample compute, once: length distribution; coverage-normalised per-position composition and Shannon entropy in 5′ and 3′ coordinates; the same partitioned per read-length class; per-position quality profile (retained at inference time even though output is FASTA — see poly-G trap, Section 7); top 3′ terminal k-mers; and a contamination k-mer screen. All later stages read from the sketch.

Stage 1 — Architecture inference (structural). Run the boundary estimators (Section 6), reconcile against any protocol hint or seqspec, and emit a candidate architecture (5′ elements, insert span, 3′ elements) with a concordance-based confidence per end and per length class.

Stage 2 — Structural self-consistency (the over/under-trim guard). After virtual trimming, assert internal coherence:

* declared constant regions are ~0 bits across their whole span (else insert was mislabelled as adapter);
* declared UMI/degenerate regions are ~2 bits (else biological sequence is about to be discarded);
* the post-trim insert terminus shows neither a residual entropy cliff (under-trim) nor loss of biological composition / a newly dominant terminal base (over-trim into the footprint);
* per-read arithmetic closes: within each length class, `read_len − adapter_onset = insert_len`, and `adapter_onset` equals the entropy-cliff position.

This stage is where trimming into the RPF body, or leaving structure on it, is caught.

Stage 3 — Identity screen (biological, provisional). On the located inserts: length-distribution shape (peaked vs broad; mode within a protocol-declared or inferred envelope, never a hardcoded gate); contamination composition by k-mer match to rRNA/tRNA/sno references; complexity/duplication, UMI-aware if a UMI was found. Emit `consistent_with_riboseq` / `inconsistent` / `indeterminate` — never "confirmed."

Stage 4 — Decision and provenance. Combine Stage 1–2 structural confidence with the Stage 3 screen into a release class; emit applied trims and the evidence object; carry `biological_confirmation: pending_alignment_gate`.

## 5. The entropy signature model

Per-position entropy in anchored coordinates is the model-free workhorse: it segments architecture without trusting an adapter database. Each element has a distinct signature:

* Constant region (adapter, linker, within-sample fixed barcode) → entropy ≈ 0 bits.
* UMI / degenerate region → ≈ 2 bits per position, base unpredictable; a high-entropy span bounded by zero-entropy constant flanks is the diagnostic UMI motif.
* Biological insert → intermediate, position-dependent (~1.6–1.9 bits), never a clean 0 or 2.

This tripartite reading turns the HMM into a principled architecture caller. States: {5′ UMI/barcode/linker, insert, 3′ adapter, 3′ UMI, poly-tail}. Emissions: per-position composition, entropy, dinucleotide. The trust gate on the HMM is: accept a segmentation only where the posterior is sharp and it agrees with the direct adapter/seqspec evidence; refuse where the posterior is diffuse.

## 6. Boundary estimators and the concordance decision function

For the 3′ insert boundary (adapter onset), per length class L, estimate the onset position by up to four independent methods:

* `b_adapter(L)` — start of an exact/fuzzy match to a known adapter (DB or seqspec). Weighted highest when present.
* `b_kmer(L)` — onset of de novo 3′ terminal k-mer enrichment (recovers adapter prefix without a DB).
* `b_entropy(L)` — the entropy-cliff position, where per-position entropy drops below a set fraction of the insert-region median.
* `b_lenmode` — onset implied by a discrete, length-correlated insert-length mode.

Reconcile into a per-class decision:

```
E        = available estimates for (end, L)
support  = reads in class L
spread   = max(E) - min(E)                      # in nt
n_indep  = number of independent estimators in E

concordance = 1.0            if spread == 0
              0.7            if spread == 1
              0.3            if spread == 2
              0.0            if spread >= 3

support_factor = min(1.0, support / S_min)      # S_min e.g. 1e4, with an
                                                # absolute floor S_floor below
                                                # which no rule may fire

confidence = concordance * support_factor
```

Decision, per end and length class:

* `safe_to_trim` — `confidence >= tau_high` AND `n_indep >= 2` AND Stage-2 self-consistency passes. For 5′ ends additionally require a named element (seqspec or high-posterior HMM segment) and `tau_high_5p > tau_high_3p`.
* `safe_to_infer_but_not_trim` — at least one estimator fires and structure is real, but concordance, support, or estimator count is below bar, or estimators disagree by ≥2 nt. No reads are modified; the sample is classified for review or protocol resolution.
* `hold` — no estimator locates a boundary, contradictory adapters/motifs coexist, or yield is near zero.

The same machinery runs for 5′ elements with the asymmetric bar. Concordance-weighted-by-support is the single scalar that becomes `trim_confidence` in the evidence object.

## 7. Checks catalogue

### Rejection checks (high specificity — where FASTQ is strong)

* Large adapter-absent fraction with a broad length distribution → not a small-RNA library. Distinguish "RNA-seq" from "read too short to reach adapter" by whether an entropy cliff exists at all.
* Flat/broad insert-length distribution after trimming → fragmentation-based, not footprints.
* Contamination k-mer screen ≈ 100% rRNA/tRNA with no plausible-length non-rRNA remainder → digestion/capture failure.
* Multiple incompatible adapters, or internal adapter k-mers → chimeras / demultiplexing failure → hold.
* No constant 3′ region and no length mode → insert cannot be located → hold, never guess.

### Artifact traps a naive extractor fails

* Poly-G dark cycles (2-colour chemistry no-calls) mimic a constant 3′ "adapter" with the same zero-entropy cliff. Discriminator: quality. A poly-G run co-occurs with a Q collapse; a real adapter does not. This is why the per-position quality profile is retained at inference time despite FASTA output. Constant 3′ + normal Q → adapter; constant 3′ + collapsed Q → base-caller artifact, flagged separately, not trimmed as adapter.
* Adapter dimers / no-insert reads — entropy cliff at position ≈ 0. Exclude; never trim to empty.
* Per-length artifacts (e.g. a non-templated 3′ base on one length class only) are invisible in pooled coordinates and require the per-length view.

### Structural self-consistency asserts (Stage 2)

The four assertions in Stage 2 above, run on virtually-trimmed reads before anything is written. Failure of any assertion demotes the sample from `safe_to_trim` to review.

## 8. Output classes

Release-compatible (all carry `biological_confirmation: pending_alignment_gate`):

* `clean_no_trim` — passes all checks; no terminal structure to remove.
* `clean_tail_depth_artifact` — low entropy only at absolute tail positions of mixed-length reads, resolved by coverage-normalised / per-length analysis; no trim applied.
* `clean_after_terminal_trim` — a strong terminal artifact was detected by convergent evidence and trimmed under explicit rules; post-trim self-consistency passes.
* `clean_after_per_length_terminal_trim` — as above, but the artifact and trim are length-specific.

Review:

* `needs_protocol_seqspec` — structured element (barcode/UMI/primer/linker/non-standard adapter) that terminal trimming cannot safely resolve.
* `needs_length_policy_review` — output dominated by short/long/disome-like fragments outside default assumptions.
* `needs_raw_fastq_review` — weak or conflicting adapter evidence, very low yield, or estimator divergence.

Hold / exclude:

* `basecaller_artifact_hold` — dominant poly-G/dark-cycle signal (quality-discriminated).
* `adapter_dimer_exclude` — insert length ≈ 0.
* `exclude_or_hold` — near-zero useful reads, unresolved contamination, or irreconcilable protocol conflict.

## 9. Evidence object

One structured record per sample (JSON + flat TSV for cohort roll-up). Extends the previous plan's schema with the boundary estimates, concordance, biological screen, quality evidence, and the provisional-confirmation flag.

```json
{
  "sample_id": "SRRXXXX",
  "getrpf_version": "x.y.z",
  "release_class": "clean_after_per_length_terminal_trim",
  "biological_confirmation": "pending_alignment_gate",
  "protocol": {
    "source": "seqspec | strong_adapter | fastqc_overrep | family_consensus | de_novo",
    "name": "Illumina TruSeq Universal",
    "adapter": "AGATCGGAAGAGCACACGTCT",
    "confidence": 0.98
  },
  "boundary_estimates_3p": {
    "by_length": {
      "29": {"b_adapter": 29, "b_kmer": 29, "b_entropy": 29, "b_lenmode": 29,
             "spread": 0, "n_indep": 4, "support": 281442,
             "concordance": 1.0, "confidence": 0.99}
    }
  },
  "trim_rules_applied": [
    {"end": "3p", "read_length_before": 29, "trim_bases": 1,
     "dominant_sequence": "A", "dominant_fraction": 0.998,
     "supporting_reads": 281442, "reason": "3p_terminal_base_bias",
     "estimators_agree": ["b_kmer", "b_entropy"]}
  ],
  "trim_rules_proposed_not_applied": [],
  "global_trim_5p": 0,
  "global_trim_3p": 0,
  "self_consistency": {"constant_region_ok": true, "umi_region_ok": null,
                       "no_residual_cliff": true, "no_overtrim": true,
                       "arithmetic_closes": true},
  "biological_screen": {
    "length_shape": "peaked",
    "insert_mode": 29,
    "contamination": {"rRNA": 0.11, "tRNA": 0.03, "other": 0.86},
    "duplication_umi_aware": 0.22,
    "verdict": "consistent_with_riboseq"
  },
  "quality_evidence": {"three_prime_q_collapse": false},
  "read_count_input": 1000000,
  "read_count_output": 812004,
  "retained_fraction": 0.812,
  "warnings": []
}
```

Cohort roll-ups: `release_qc_summary.tsv`, `release_qc_flags.tsv`, `adapter_protocol_summary.tsv`, `samples_needing_seqspec.tsv`, `samples_excluded_or_held.tsv`.

## 10. Execution modes and entry point

Production — infer on subsample, apply frozen rules on the stream, emit cleaned RPF FASTQ + collapsed FASTA + evidence object + applied/proposed rules + release class. High-confidence safe trims are applied automatically; every edit is written with the evidence that triggered it. No edit is hidden.

Audit (`--audit-only`) — same inference, no application. Emits proposed rules and the diagnostics explaining which would have fired. For rule development and new-family validation; not the portal path.

Override — apply explicit user-supplied rules; emit cleaned output, provenance, and before/after QC. For reproducing a frozen release or handling a known edge case.

Samplesheet entry point — accept local FASTQ paths so extraction runs without re-downloading or fishing inputs out of old work directories.

```
sample_id,fastq_1,fastq_2,layout,organism,study_id,protocol_hint,fastqc_dir   # full
sample_id,fastq_1                                                              # minimum
```

Reuses the same modules as the download path; stops after release-ready collapsed FASTA + provenance.

Per-sample outputs: `{s}.rpfs.fastq.gz`, `{s}.rpfs.collapsed.fa.gz`, `{s}.evidence.json`, `{s}.terminal_diagnostics.tsv`, `{s}.per_length_terminal_diagnostics.tsv`, `{s}.trim_rules.applied.json`, `{s}.trim_rules.proposed.json`, `{s}.release_qc.json`.

## 11. Command-line surface

The refined CLI is verb-first and maps onto the algorithm rather than onto historical accretion. The two facts it has to make ergonomic are that inference runs on a subsample while application streams, and that refusal is a normal outcome that a pipeline must be able to branch on without parsing JSON.

```
getRPF
  extract   INPUT -o DIR          # production path: sketch → infer → self-consistency → apply → screen → provenance
      --protocol-hint NAME         #   or  --seqspec-dir DIR  |  --architecture-db JSON
      --fastqc-dir DIR             # optional FastQC evidence to fold into protocol assignment
      --infer-reads N              # subsample depth for Stage 0/1 (default 5e5); apply still streams all reads
      --audit-only                 # infer + screen, apply nothing (rule development / new-family validation)
      --rules FILE                 # override mode: apply explicit frozen rules instead of inferred ones
      --fail-on {none,review,hold} # process exit-code threshold (default: none)
      --collapse / --no-collapse
  sketch    INPUT -o FILE          # Stage 0 only, reference-free: length dist; per-position & per-length
                                   #   5′/3′ entropy + composition; quality profile; 3′ terminal k-mers;
                                   #   contamination k-mer screen. The substrate everything else reads from.
  run       SAMPLESHEET -o DIR     # cohort driver over local FASTQ; per-sample outputs + cohort roll-ups
                                   #   (accepts the same inference/apply flags as extract)
  seqspec
      generate  INPUT              # de novo architecture → candidate seqspec
      promote   --family GLOB      # promote a recurring de novo structure to a named architecture
      validate  SEQSPEC INPUT      # check a seqspec against a sample's sketch
  ingest    EVIDENCE... --db qc.db # provenance JSON/TSV → DuckDB review store
  plot
      entropy   INPUT              # per-position 5′/3′ entropy with HMM segments
      terminal  INPUT              # per-length terminal diagnostics
      softclips --align-json J     # alignment soft-clip summary (second gate)
  align-detect INPUT --star-index IDX   # SECOND GATE: biological confirmation. Not a trimming input.
```

Modes are flags, not verbs. Production, audit, and override are the same inference path with application toggled (`--audit-only`, `--rules`). Keeping them as one command is what lets the M6 acceptance criterion — audit and production agree on proposed rules — hold by construction rather than by luck.

Refusal is machine-readable at the process level. `extract` and `run` always write the release class to `release_qc.json`, and additionally mirror it to the exit code so a Nextflow or CI step can branch without reading files: `0` = all processed samples release-compatible, `10` = at least one review class, `20` = at least one hold/exclude. `--fail-on` sets where nonzero begins (e.g. `--fail-on hold` treats review as acceptable). The exit code is a convenience mirror, never the sole record.

The FASTQ stage stays reference-free. `extract` never touches an aligner; the biological gate is a separate `align-detect` invocation by design, so this stage runs anywhere and stays independently auditable. A `--gate align` convenience may chain the two, but decoupled is the default.

Consolidation from the current surface. `decide-trim` folds into `extract` — its architecture/alignment consensus becomes one estimator feeding the concordance function of Section 6, not a standalone decision. `extract-rpf` becomes an alias of `extract`; `check`/`check-cleanliness` are replaced by `sketch` plus the release classifier, since the pass/fail cleanliness report is superseded by the evidence object; `plot-hmm`→`plot entropy`, `plot-softclips`→`plot softclips`, `ingest-duckdb`→`ingest`. `detect-adapter` is retained as a targeted single-adapter probe but is no longer the basis of any trimming decision. Deprecated spellings survive one release as hidden aliases emitting the new outputs.

## 12. Handoff to the alignment gate

Reserved for the alignment stage, and deliberately not attempted here: triplet periodicity, CDS enrichment, start/stop metagene shape, P-site-offset-by-length, soft-clip asymmetry, and RiboMetric-style biological QC. These earn the biological claim.

One tempting boundary case: a pseudo-periodicity signal is obtainable by pseudo-aligning inserts to a frame-annotated transcriptome k-mer set without full alignment. Resist it. It smuggles the alignment gate forward, couples this stage to a reference, and blurs the structural/biological split that makes each stage independently auditable.

The sharpest pre-alignment biological screening axis, and the one to lead with, is the joint distribution of insert length and adapter-onset discreteness: footprint libraries have inserts that are short, tightly length-distributed, and almost all carry a 3′ adapter beginning at a discrete, length-correlated position. RNA-seq fails the discreteness; random degradation fails the tightness. It does not separate Ribo-seq from other tight small-RNA protocols — nothing pre-alignment does — but it falls straight out of the Stage-0 sketch.

## 13. Milestones

M0 — Entropy consolidation. Retire `check.py`'s bespoke entropy; route cleanliness QC through `SignalProcessor`. Acceptance: a max-complexity, low-coverage tail position no longer fails; mixed 28–32 nt clean RPFs pass; synthetic 3′ A=99% reads flag as `three_prime_terminal_bias`. Smallest change, removes the motivating false positives — ship alone.

M1 — Sketch + coordinate/per-length diagnostics. Extend `SignalProcessor` to emit the Stage-0 sketch including per-length 5′/3′ terminal diagnostics and the quality profile. Acceptance: per-length terminal artifact is visible where pooled analysis misses it; poly-G run is distinguishable from adapter by Q.

M2 — Boundary estimators + concordance + refusal. Implement the four 3′ estimators (and 5′ variants), the concordance/confidence function, and the three-way decision with the asymmetric 5′ bar. Acceptance: synthetic 3′ 1-base artifact yields exactly one 3′ trim; synthetic internal motif proposes nothing; deliberately divergent estimators route to `needs_raw_fastq_review` rather than trimming.

M3 — Apply path + provenance. Thread per-length override trims into `RPFExtractor.extract_rpfs`; connect the M2 decision (and the existing `TrimDecider` consensus) to actual application; emit `trim_rules.applied.json`. This is the real unlock. Acceptance: applied rules are exactly reproducible from provenance; before/after shows the terminal artifact removed and the insert body untouched.

M4 — Identity screen. Length-shape classifier, rRNA/tRNA/sno k-mer contamination screen, UMI-aware duplication, poly-G/dimer traps. Acceptance: an RNA-seq library and a degradome library are both screened out with the correct distinguishing reason; a clean footprint library returns `consistent_with_riboseq`.

M5 — Release classifier + evidence object + cohort store. Emit the Section-8 classes and Section-9 evidence object; rewrite `categorize_failures` so it preserves subclasses rather than re-collapsing them; roll up into DuckDB and the cohort TSVs. Acceptance: every released FASTA has an `evidence.json`; every held sample has a machine-readable reason; the cohort summary is reproducible from emitted TSV/JSON.

M6 — Samplesheet + audit/override modes. Local-FASTQ driver over existing modules; `--audit-only` and override flags. Acceptance: a local-path samplesheet processes a subset with no ENA/SRA download; audit and production agree on proposed rules.

M7 — Portal release audit. Run the 6k local subset and 4k cluster set; stratify by protocol and class; run the alignment gate on representative and all suspicious groups; freeze release classes and exclusions. Acceptance: release/exclusion decisions are reproducible from provenance; the seqspec promotion queue is populated from recurring de novo families.

Sequencing note: breadth of `safe_to_trim` comes chiefly from protocol assignment (Rule-B/D territory), not from the terminal-bias engine, so the seqspec promotion loop runs alongside M4–M5 and is where recall actually grows once M0–M3 have made the structural decisions safe.

## 14. Non-goals and limits

* No biological confirmation of Ribo-seq identity from FASTQ alone. The strongest positive verdict is provisional.
* No reference genome/transcriptome alignment, and no pseudo-alignment periodicity, at this stage.
* No hardcoded RPF length window as a gate; length envelopes are protocol-declared or inferred.
* No trimming on single-line evidence, and no 5′ trimming on de novo terminal bias alone.
* Manual seqspec authoring is a fallback, not the workflow; the intended path is automatic protocol assignment backed by a growing, promotion-fed architecture library.
