# Real-library panel v1: P6 audit (2026-09-13)

- **Code:** branch `feat/read-structure-inference` at `a71bc3d`, transform schema `getrpf.structure.transform/3`.
- **Data:** the first 1M reads of each run (`scripts/fetch_panel_subsets.sh`). Inference uses reads 1–300k.
- **Command:** `scripts/run_structure_panel.py --data-dir validation/panel_v1/data --output-dir validation/panel_v1/results` (audit only). The files in this directory are its output; this README was added afterwards.
- **No production FASTQ was written or promoted.**

## Result

| Run | Architecture | Transform | Audit, reads 1–300k | Held out, reads 600k–1M |
|---|---|---|---|---|
| SRR1944950 | `[nta 0-1 kept][insert]` | emitted | 299,987 accepted (100.0%) | 100.0% accepted; mode 31 nt; 99.8% at 26–34 nt; Q5 `riboseq_likely` |
| SRR3945920 | `[insert][nta 0-1 kept][adapter AGATCGGAAGAG...]` | emitted | 293,296 (97.8%) | 97.7%; mode 30 nt; 95.1%; `riboseq_likely` |
| SRR3945930 | `[random 3 UMI][nta 0-1 kept][insert][random 4 UMI][adapter CTGTAGGCACCA...]` | emitted | 290,420 (96.8%) | 96.8%; mode 29 nt; 93.0%; `riboseq_likely` |
| SRR12693498 | `[random 13 UMI][fixed GGG][insert][nta 0-1 kept][poly(A)][adapter AGATCGGAAGAG...]` | emitted | 256,443 (85.5%) | 85.5%; mode 33 nt; 74.5%; `riboseq_likely` |
| SRR23242345 | `[random 5 UMI][nta 0-2 kept][insert][adapter AGATCGGAAGAG...]` | **withheld** | — | — |

SRR23242345 is withheld because read positions 6–7 look non-templated in about 78% of reads (spec §7.2).

## What counts as evidence

**Independent.** These expectations were fixed before analysis in `tests/test_structure/test_real_panel.py` and `test_real_extraction.py`:
- 4 of 5 expected structures confirmed.
- The OTTR claim of a 7 nt UMI is contradicted: positions 6–7 are skewed (60% T, 78% C), not random. It is kept as a strict expected failure.
- All 5 architectures are reproduced from held-out reads 300k–600k.
- Q5 is `riboseq_likely` on held-out reads 600k–1M for all 4 emitted runs.
- Audit and production agree on real reads.

**Regression guard only.** `benchmark.md` passes (layout 5/5, adapter 5/5, tail 5/5, transform 5/5, UMI 3/3, RPF recovery 4/4, no false-safe transforms). Its expectations in `truth.yaml` were written from the same measurements, so a pass shows the method still reproduces them; it is not independent validation.

## Audit rejections (reads 1–300k)

| Run | Anchor not found | Adapter dimer | Fixed-block mismatch | Shorter than 20 nt | Longer than 40 nt | No insert |
|---|---|---|---|---|---|---|
| SRR1944950 | – | – | – | – | 13 | – |
| SRR3945920 | 3,185 | 367 | – | 3,135 | 17 | – |
| SRR3945930 | 6,065 | 980 | – | 2,535 | – | – |
| SRR12693498 | 29,038 | 3 | 1,934 | 9,239 | 3,047 | 296 |

SRR12693498 boundary counts: 268,996 at the poly(A) tail and 29 at the adapter. They reconcile exactly: 256,443 accepted + 12,582 rejected by length or empty insert = 269,025.

The earlier uncommitted decision record, in the git-ignored `validation_rewire_v1_hpc_retry/`, gave 89,634 + 19 and did not reconcile. That record is superseded by this README.

## Known limitations and open items

1. **Depth fragility, SRR3945930.** At 50k and 100k reads, read position 3 has within-fragment agreement 0.43 and 0.41, just above `agree_random` (0.40). Q2 then becomes `[random 2][nta 0-2]` and the transform is withheld: it fails closed, but it is not stable. At 200k and 300k reads the value is 0.38 and 0.39, clearing the threshold by 0.01. The random-level test needs a depth-aware or relative criterion before production use.
2. **SRR1944950 5′ NTA is unconfirmed.** Confirming or refuting it is phase P7, by alignment. NTA rate estimates are not calibrated.
3. **Q5 is length-only.** Frame periodicity needs alignment.
4. **Leading reads only.** The subsets come from the first flow-cell tiles.
5. **SRR23242345 (OTTR).** Re-read the methods for read positions 6–7 before deciding whether to remove them.
6. **Missing panel case.** A McGlincy-style run (3′ UMI plus inline barcode) is still to be sourced.

## Default extraction path (P6 decision)

**Recommendation:** keep the legacy `extract` behaviour as the default. The new path is available now as `getRPF infer-structure` followed by `getRPF extract --architecture <sample>.seqspec.yaml`. Make it the default after:
- P7 has confirmed or refuted the SRR1944950 NTA;
- item 1 is fixed;
- a full-depth production run on HPC shows audit and production agreeing over whole files;
- the McGlincy-style case is in the panel.

This decision belongs to the maintainer.
