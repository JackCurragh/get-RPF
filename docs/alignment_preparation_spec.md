# getRPF trimming and alignment preparation

**Status** Draft for implementation, 23 November 2026

Produce enough correctly positioned ribosome footprints for reproducible profiles. getRPF applies terminal cuts supported by the sequence evidence, then sends the resulting reads to the main alignment. The pipeline uses annotated CDS phase, read depth, and length-specific offsets to decide which alignments contribute to profiles.

This specification covers single-end Ribo-seq in the per-sample production workflow. It builds on the [v1 read-structure model](read_structure_inference_spec.md). The alignment input carries candidate footprints. A complete molecular architecture remains a stronger claim that getRPF can report when supported.

## Evidence behind the changes

| Finding | Implementation consequence |
| --- | --- |
| For `SRR836195`, the first million reads supplied 39 informative groups. A million reads dispersed across the same run supplied 1,081 and resolved the trim. | Retry unresolved inference on a representative sample. |
| In `SRR392986`, 47,624 of 48,383 held-out reads with the reported adapter motif carried `ATC` immediately upstream. getRPF found the longer motif. Cutting at the shorter motif produced a modal 3-nt soft clip. Q3 stayed unresolved when the search window increased from 16 to 32 nt. | Verify the insert-facing edge of supplied motifs and apply a supported adapter cut even when adjacent bases remain uncertain. |
| At the minimum scoring depth of three other sequences, full credit for a tied consensus produces 43.75% expected agreement for uniform random bases. The random cutoff is 40%. A separate synthetic probe showed that an intermediate-agreement base can be cut when a deeper position is classified as fixed. | Fix tie handling and require evidence for the full inferred trim span. |
| On 50,000 held-out reads, the reported trim preserved the strong `SRR1062396` profile and recovered 25,607 genome placements for `SRR392986`. The latter had 37% dominant CDS frame share and zero recommended lengths. | Evaluate recovery using supported footprint lengths and offsets as well as alignment yield. |

The local pilot used the frozen zebrafish reference and the same 50,000 held-out reads per recipe. The workspace has the underlying read traces in `tmp/local-read-trace.jodJ5b` and the reported recipe comparison in `tmp/reported-ribo-trims.4mr3zn2o/comparison.tsv`.

## Run recipes

Add an optional Nextflow parameter `getrpf_trim_recipes`. Read the CSV once and join it by raw `run_id` before library aggregation. Accept the supplied header `run_id,adapter,5p_barcode_length,3p_barcode_length`. Support optional `source` and `input_state` columns for curated entries. Set `input_state` to `pretrimmed` only when the deposited input state has been verified. Trim lengths are non-negative integers and zero is an explicit value. A blank field is unspecified. Treat `passed` as a status marker with unspecified adapter sequence unless its meaning has been verified. Report malformed values and conflicting duplicate accessions with the source row.

| Evidence | Action |
| --- | --- |
| A supplied adapter and proposed trim agree with the observed layout | Use the supplied values for the supported parts of the read. Record their source in the extraction report. |
| The reads support a motif extending upstream of the supplied adapter sequence | Use the supported cut position and record the difference from the supplied recipe. |
| A field is missing or contradicted by the reads | Infer that part of the structure. Record the reason when a reported cut is declined. |
| Provenance and read evidence support a pre-trimmed input | Use the deposited end for the corresponding boundary. |

Verify each action on the inference sample using the existing anchor search and terminal profile. Check the orientation, insert-facing adapter start, and full extent of each proposed cut. A read-length mode can guide the search. A fixed trim length needs sequence evidence covering that span. Use the same supported cut rules for curated and inferred blocks.

## Inference changes

1. **Retry only when a pending decision could change a cut.** Start with the current prefix sample. If an adapter or junction stays unresolved, retry once using a deterministic, bounded sample drawn across the whole FASTQ. Use a streaming reservoir for sequential gzip input. Record sample size, sampling method, and seed. Where two adequately supported calls disagree, retain the contested span and report both calls.
2. **Retain evidence independently for both termini.** Keep adapter location, the 5′ terminal block, and the 3′ pre-adapter block as separate decisions. Build a 5′ profile from reads with a credible grouping frame even when the 3′ endpoint is unavailable. For confidently located anchors, use the existing adapter matcher when streaming the full FASTQ. Short terminal adapter matches require a supported motif and a per-read end match.
3. **Fix consensus support.** Give a base matching one of `k` tied consensus bases `1/k` agreement credit. Keep distinct-sequence depth and scored position counts in the report. Use the existing benchmark to check calibration and recovery before changing the group minimum or agreement cutoffs.
4. **Apply only supported terminal spans.** Start at the relevant edge and stop at the first position without sufficient evidence for removal. A supported compound motif can establish the whole span. Preserve ambiguous junction bases, and preserve the existing cut convention for supported tails. Retain per-position evidence when a cut is ambiguous or underpowered.

Implement these changes in `core/structure/observe.py`, `anchors.py`, `pileup.py`, and `assemble.py`. Keep the adapter and each junction answer available to the stream writer even when assembly cannot produce a complete architecture.

## getRPF command and output contract

Add `getRPF extract --prepare-for-alignment` to the structure extraction path. An optional `--trim-recipe <json>` supplies one normalized row with the same field names as the CSV, integer trim lengths, and nulls for unspecified fields. Reuse the adapter locator, fixed-block checks, sequence and quality slicing, and FASTA collapse for full and partial cuts. Make the fragment length policy a separate step from slicing so preparation can defer it until after alignment. Share the report and collapsed-FASTA writer with the explicit `--architecture` path when preparation is requested. Existing strict callers can continue to request the v1 transform contract. The entry points are `cli.py` and `core/structure/workflow.py`. Keep the per-read operation in `core/structure/transform.py`.

For each valid record, resolve a zero-based, half-open interval `[start, end)` on the original read. Start with the complete input sequence. Advance `start` for a supported 5′ block whose per-read checks pass. Set `end` at a supported adapter match, then move it inward for a supported adjacent block. A verified pre-trimmed read end can also anchor a 3′ block. If a cut lacks evidence for that read, retain those bases. Emit unresolved records with supported cuts, or with the full input sequence when neither cut is supported. Record positive adapter dimers and empty intervals as read rejections.

Emit the existing collapsed FASTA with one entry per sequence and an `_xN` multiplicity. Keep the report’s structure confidence under `transform`. Add a `preparation` section with counts of emitted reads by endpoint status, the count of unmodified reads, endpoint evidence, actual bases cut per end, unmatched-adapter counts, and recipe provenance. Use the existing extraction section for lengths and rejections and the existing reads field for the output name. Set `preparation.requires_local_alignment` when any emitted read has an unresolved end. A library with unresolved ends and emitted reads has successful alignment preparation.

Every input read must be accounted for. On a successful run, the sum of the FASTA multiplicities plus the sum of recorded read rejections equals the input count. A complete FASTQ parse is required for publication. Preserve sequence and quality alignment when writing FASTQ. Treat random-block removal and UMI deduplication as distinct operations. The collapsed alignment input represents raw read multiplicity.

## Nextflow integration

1. Extend `EXTRACT_RPFS` to receive the matched recipe, invoke preparation, and publish the same required extraction report and collapsed FASTA channels. Pass the selected alignment mode downstream as sample metadata with a stable join key.
2. Route prepared reads through the existing per-sample alignment workflow. Use Local STAR when the output contains unresolved ends and use the configured alignment mode for fully resolved output. Run STAR once per run and derive the genome and transcriptome views from that invocation. In preparation mode, compute any clipping evidence from this BAM in place of the extraction-time STAR check.
3. Defer the 20–40 nt footprint filter for prepared input to the aligned-footprint eligibility step. The query length before alignment can include technical bases. Publish counts for input, mapped, and eligible reads with multiplicities restored. Treat pre-alignment cleanliness findings as diagnostic metrics.
4. Retire the raw hold FASTQ, `ALIGNMENT_GETRPF_HOLDS`, its STAR process and publishers, and `align_getrpf_holds` after the replacement path passes integration tests. Use the production BAM and extraction report for troubleshooting.

The integration starts in `pipelines/riboseq/subworkflows/data_acquisition.nf`, `modules/extract_rpfs.nf`, `subworkflows/alignment.nf`, and `modules/star_align.nf` in `ensembl-genes-nf`. Extend the matrix workflow after validating the per-sample contract and its count provenance.

## Footprint coordinates and profile eligibility

Use the same CIGAR interpretation in RiboMetric offset estimation, release QC, and BAM-to-track projection. Today RiboMetric can subtract leading soft clips from the aligned start, while track generation indexes from the first aligned reference base. A clipped technical prefix makes those two coordinate origins disagree. The implementation points are `RiboMetric/bam_processing.py` and `RiboMetric/qc.py` in the RiboMetric checkout and `pipelines/riboseq/bin/bam_to_bed.py` in the pipeline checkout.

Define the aligned query length as the sum of `M`, `=`, `X`, and `I`. Count the 5′-relative offset from the first aligned query base in read orientation. Walk CIGAR to project that offset to the reference. Insertions consume a query offset and have no reference coordinate at that offset. Deletions and introns advance the reference. Hard and soft clips contribute zero query bases to the footprint. Apply the same eligibility predicate before estimating length-specific offsets and before generating tracks. Test offset projection in both orientations and across splice junctions.

For the first release, treat any remaining 5′ soft clip as an unresolved biological origin. Keep those alignments in the main BAM and mapped counts. Offset estimation and profile generation require an aligned 5′ end. A 3′ clip can remain outside the aligned span if the resulting length and 5′ origin satisfy footprint eligibility. Apply the 20–40 nt policy to the aligned query length and account explicitly for any alignment the projection cannot place. Reuse the same predicate with alternate offset providers before enabling them for prepared reads.

Use existing production depth, supported-length and offset gates. Report both the share of mapped reads eligible for profiling and the share assigned to supported lengths. Evaluate candidate recipes on the same held-out inputs. Estimate offsets on one set and verify expected CDS frame and start/stop profiles on independent reads or genes. Compare supported depth together with the phase distribution and offset stability.

## Build order and acceptance

| Step | Scope | Acceptance |
| --- | --- | --- |
| 1 | getRPF inference and CLI | Exercise provided, inferred, zero-trim, and unresolved cases. The uncertain-base and tied-consensus probes become regression tests. Reads accepted by the strict transform in fully supported controls reproduce their existing trimmed sequence and multiplicities. Reconcile additional prepared reads separately. |
| 2 | Alignment footprint geometry | Use shared fixtures for `28M`, `3S28M`, `28M4S`, reverse-strand cases, an insertion, a deletion, and a splice. RiboMetric and track generation agree on eligible lengths, offset origins and projected sites. Include lengths 20 and 24 to catch any legacy consumer minimum. |
| 3 | Per-sample Nextflow wiring | Run fixtures for a matched recipe, absent recipe, unsupported recipe, and unresolved end. Confirm one main STAR invocation, required reports, coherent counts, and publication through normal QC. Confirm that removing the diagnostic branch leaves the main execution graph valid. |
| 4 | Paired zebrafish pilot | Freeze the FASTQs, reference, validation partitions and QC rules. Compare current extraction with preparation on `SRR1062396` and `SRR13998266` as passing controls, `SRR836195` as the sampling case, `SRR392986` as the adapter case, and `SRR21218374` as a low-periodicity control. Expand to the remaining held cohort after these checks pass. |

The pilot succeeds when passing controls retain their supported profiles and recovery produces additional reads at demonstrably supported lengths. Report any run that gains alignments without gaining supported footprint depth as a mapping gain. Keep the existing release gate in force for that run.

Run the structure suite from the get-RPF checkout.

```bash
uv run --no-project --python 3.12 --with . --with pytest \
  pytest -o addopts='' -p no:cacheprovider -q \
  tests/test_structure/test_{anchors,pileup,assemble,workflow,transform,align}.py
```

Run the pipeline QC regression tests from the ensembl-genes-nf checkout.

```bash
uv run --no-project --python 3.12 --with pytest --with pyyaml \
  pytest -q -p no:cacheprovider pipelines/riboseq/tests/test_qc_gate_metrics.py
```
