"""Step 4: reference alignment, as an escalation and as a check (spec §5.4).

One local STAR alignment of a bounded set of sequences. It records, in read
orientation, the 5' and 3' soft-clip length histograms, the clipped bases and
the mismatch rate inside the alignment. It never compares mapping rates
across candidate trims.

Aligning the inserts an architecture would emit tests the pileup's Q2/Q3
answers against the genome:
- a technical block the pileup missed shows as clips of 2 nt or more;
- a real non-templated junction base shows as a 1 nt clip spike, well above
  the error rate, whose bases match the pileup's disagreeing bases.
Alignment never silently overrides the pileup. A missed technical block makes
the answer ``conflicting``: the structure itself is in doubt, so the transform
is withheld. A disagreement about one junction base only flags the transform,
because that base is kept in the insert either way (spec §7.2).
"""

from __future__ import annotations

import gzip
import math
import re
import shutil
import subprocess
import tempfile
from collections import Counter
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .assemble import StructureInference, build_architecture, decide_transform
from .config import InferenceConfig
from .model import Answer, Architecture, Evidence, JunctionCall, Status
from .transform import Accepted, Transform

_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
_MD = re.compile(r"(\d+)|(\^[A-Z]+)|([A-Z])")
_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


@dataclass
class AlignmentProfile:
    sequences: int = 0
    aligned: int = 0
    clip5: Counter[int] = field(default_factory=Counter)
    clip3: Counter[int] = field(default_factory=Counter)
    clipped5_bases: Counter[str] = field(default_factory=Counter)
    """Outermost clipped base at the 5' end, in read orientation."""
    clipped3_bases: Counter[str] = field(default_factory=Counter)
    aligned_bases: int = 0
    mismatches: int = 0

    @property
    def mismatch_rate(self) -> float:
        return self.mismatches / self.aligned_bases if self.aligned_bases else 0.0

    def fraction(self, end: str, low: int, high: Optional[int] = None) -> float:
        clips = self.clip5 if end == "5" else self.clip3
        hits = sum(
            count
            for length, count in clips.items()
            if length >= low and (high is None or length <= high)
        )
        return hits / self.aligned if self.aligned else 0.0


@dataclass(frozen=True)
class JunctionCheck:
    question: str
    verdict: str
    """confirmed | conflicting | inconclusive | clean"""
    one_base_clip: float
    longer_clip: float
    longer_clip_mode: Optional[int]
    clipped_bases: Dict[str, float]
    answer: Answer
    aligned: int = 0
    """Unique alignments behind the verdict."""
    flag: Optional[str] = None
    """Set when the alignment disagrees about a junction base: reported on the
    transform, never withholding it."""


def parse_sam(lines: Iterable[str]) -> AlignmentProfile:
    """Clip and mismatch statistics from primary alignments in SAM text."""
    profile = AlignmentProfile()
    for line in lines:
        if line.startswith("@"):
            continue
        fields = line.rstrip("\n").split("\t")
        flag = int(fields[1])
        if flag & (0x4 | 0x100 | 0x800):
            continue
        ops = [(int(n), op) for n, op in _CIGAR.findall(fields[5])]
        if not ops:
            continue
        left = ops[0][0] if ops[0][1] == "S" else 0
        right = ops[-1][0] if ops[-1][1] == "S" else 0
        if flag & 0x10:
            clip5, clip3 = right, left
            read = fields[9][::-1].translate(_COMPLEMENT)
        else:
            clip5, clip3 = left, right
            read = fields[9]
        profile.aligned += 1
        profile.clip5[clip5] += 1
        profile.clip3[clip3] += 1
        if clip5:
            profile.clipped5_bases[read[0]] += 1
        if clip3:
            profile.clipped3_bases[read[-1]] += 1
        profile.aligned_bases += sum(n for n, op in ops if op in "M=X")
        md = next((tag[5:] for tag in fields[11:] if tag.startswith("MD:Z:")), "")
        profile.mismatches += sum(1 for _, _, base in _MD.findall(md) if base)
    return profile


def run_star(
    sequences: Sequence[str],
    star_index: Path,
    config: Optional[InferenceConfig] = None,
) -> AlignmentProfile:
    """Align sequences locally with STAR (unique alignments only)."""
    config = config or InferenceConfig()
    if not sequences:
        return AlignmentProfile()
    star = shutil.which("STAR")
    if star is None:
        raise RuntimeError("STAR not found on PATH")
    with tempfile.TemporaryDirectory(prefix="getrpf_align_") as directory:
        workdir = Path(directory)
        fastq = workdir / "reads.fastq"
        with fastq.open("w") as handle:
            for index, sequence in enumerate(sequences):
                handle.write(f"@s{index}\n{sequence}\n+\n{'I' * len(sequence)}\n")
        subprocess.run(
            [
                star,
                "--runMode", "alignReads",
                "--genomeDir", str(star_index),
                "--readFilesIn", str(fastq),
                "--outFileNamePrefix", f"{workdir}/",
                "--outSAMtype", "SAM",
                "--outSAMattributes", "NH", "MD",
                "--outSAMunmapped", "None",
                "--alignEndsType", "Local",
                "--outFilterMultimapNmax", "1",
                "--outFilterMismatchNmax", "3",
                "--outFilterScoreMinOverLread", "0",
                "--outFilterMatchNminOverLread", "0",
                "--outFilterMatchNmin", str(config.align_min_matched),
                "--runThreadN", str(config.align_threads),
            ],
            check=True,
            capture_output=True,
            text=True,
        )  # fmt: skip
        if _star_input_reads(workdir / "Log.final.out") == 0:
            # Seen with the bioconda osx-arm64 STAR 2.7.11b build: it exits 0
            # having read nothing, which would otherwise read as "nothing
            # aligns".
            raise RuntimeError(
                f"STAR read 0 of {len(sequences)} input reads; this STAR build "
                "cannot read FASTQ input"
            )
        with (workdir / "Aligned.out.sam").open() as sam:
            profile = parse_sam(sam)
    profile.sequences = len(sequences)
    return profile


def _star_input_reads(log: Path) -> Optional[int]:
    """'Number of input reads' from STAR's Log.final.out, when present."""
    if not log.exists():
        return None
    for line in log.read_text().splitlines():
        if "Number of input reads" in line:
            return int(line.split("|")[1])
    return None


def build_star_index(fasta: Path, out_dir: Path, threads: int = 4) -> Path:
    """Build a STAR index for a small genome (yeast-sized)."""
    star = shutil.which("STAR")
    if star is None:
        raise RuntimeError("STAR not found on PATH")
    out_dir.mkdir(parents=True, exist_ok=True)
    genome = out_dir / "genome.fa"
    opener = gzip.open if str(fasta).endswith(".gz") else open
    length = 0
    with opener(fasta, "rt") as source, genome.open("w") as target:  # type: ignore[operator]
        for line in source:
            target.write(line)
            if not line.startswith(">"):
                length += len(line.strip())
    sa_bases = min(14, int(math.log2(max(length, 2)) / 2 - 1))
    subprocess.run(
        [
            star,
            "--runMode", "genomeGenerate",
            "--genomeDir", str(out_dir),
            "--genomeFastaFiles", str(genome),
            "--genomeSAindexNbases", str(sa_bases),
            "--runThreadN", str(threads),
            "--outFileNamePrefix", f"{out_dir}/",
        ],
        check=True,
        capture_output=True,
        text=True,
    )  # fmt: skip
    return out_dir


def emitted_inserts(
    architecture: Architecture,
    reads: Sequence[str],
    config: Optional[InferenceConfig] = None,
) -> List[str]:
    """The inserts the transform would emit, junction bases kept (v1)."""
    config = config or InferenceConfig()
    transform = Transform(architecture, config)
    inserts: List[str] = []
    for read in reads:
        outcome = transform.apply(read, "I" * len(read))
        if isinstance(outcome, Accepted):
            inserts.append(outcome.insert)
            if len(inserts) == config.align_reads:
                break
    return inserts


def check_junction(
    answer: Answer, end: str, profile: AlignmentProfile, config: InferenceConfig
) -> JunctionCheck:
    """Compare one pileup junction answer (Q2: end '5', Q3: end '3') with
    the alignment of the emitted inserts."""
    one = profile.fraction(end, 1, 1)
    longer = profile.fraction(end, 2)
    clips = profile.clip5 if end == "5" else profile.clip3
    longer_mode = max(
        (length for length in clips if length >= 2),
        key=lambda length: clips[length],
        default=None,
    )
    bases = profile.clipped5_bases if end == "5" else profile.clipped3_bases
    total = sum(bases.values())
    composition = {b: round(bases[b] / total, 3) for b in "ACGT"} if total else {}
    top_clipped = (
        max(composition, key=lambda b: composition[b]) if composition else None
    )
    call = answer.value if isinstance(answer.value, JunctionCall) else None
    pileup_top = (
        max(call.nta_bases, key=lambda b: call.nta_bases[b])
        if call is not None and call.nta_bases
        else None
    )
    error = profile.mismatch_rate
    side = "5'" if end == "5" else "3'"

    if profile.aligned < config.align_min_aligned:
        # Too few alignments is missing evidence, never a contradiction.
        verdict = "underpowered"
        text = (
            f"only {profile.aligned} of {profile.sequences} inserts aligned "
            f"uniquely (need {config.align_min_aligned}); no verdict"
        )
    elif longer >= config.align_missed_block_rate:
        verdict = "conflicting"
        text = (
            f"alignment shows {side} clips of >=2 nt in {longer:.0%} of aligned "
            f"inserts (mode {longer_mode} nt): a technical block the pileup did "
            "not call"
        )
    elif call is not None and call.nta_length:
        if (
            one >= config.align_confirm_rate
            and one >= config.align_error_multiple * error
            and (pileup_top is None or top_clipped == pileup_top)
        ):
            verdict = "confirmed"
            text = (
                f"alignment confirms a {side} non-templated base: 1 nt clips in "
                f"{one:.0%} of aligned inserts (mismatch rate inside alignments "
                f"{error:.2%}), clipped bases {composition.get(top_clipped or '', 0):.0%} "
                f"{top_clipped}, as the pileup found"
            )
        elif one < config.align_refute_rate:
            verdict = "conflicting"
            text = (
                f"alignment does not support the {side} junction base: 1 nt clips "
                f"in only {one:.1%} of aligned inserts"
            )
        else:
            verdict = "inconclusive"
            text = (
                f"alignment is inconclusive about the {side} junction base: 1 nt "
                f"clips in {one:.0%} of aligned inserts, clipped bases mostly "
                f"{top_clipped} (pileup: {pileup_top})"
            )
    elif one >= config.align_confirm_rate:
        verdict = "conflicting"
        text = (
            f"alignment finds {side} 1 nt clips in {one:.0%} of aligned inserts, "
            "where the pileup called no junction base"
        )
    else:
        verdict = "clean"
        text = (
            f"alignment agrees: {side} clips of 1 nt in {one:.1%} and >=2 nt in "
            f"{longer:.1%} of aligned inserts"
        )

    evidence = answer.evidence + (
        Evidence(
            "alignment", f"clip{end}_1nt_fraction", round(one, 4), profile.aligned
        ),
        Evidence(
            "alignment",
            f"clip{end}_2plus_fraction",
            round(longer, 4),
            profile.aligned,
            note=f"mode {longer_mode} nt" if longer_mode else "",
        ),
        Evidence("alignment", f"clipped{end}_bases", composition, total),
        Evidence(
            "alignment",
            "mismatch_rate",
            round(error, 5),
            profile.aligned_bases,
            note="mismatches per aligned base inside local alignments",
        ),
        Evidence("alignment", "verdict", verdict, profile.aligned, note=text),
    )
    # Only a missed technical block puts the structure in doubt. A junction
    # base stays in the insert whichever method is right, so a disagreement
    # about it flags the transform instead of withholding it (spec §7.2).
    blocking = verdict == "conflicting" and longer >= config.align_missed_block_rate
    status = Status.CONFLICTING if blocking else answer.status
    flag = (
        f"{answer.question}: {text}"
        if verdict == "conflicting" and not blocking
        else None
    )
    updated = replace(
        answer,
        evidence=evidence,
        status=status,
        explanation=f"{answer.explanation} Alignment: {text}.",
    )
    return JunctionCheck(
        answer.question,
        verdict,
        one,
        longer,
        longer_mode,
        composition,
        updated,
        profile.aligned,
        flag,
    )


def apply_alignment_check(
    result: StructureInference,
    reads: Sequence[str],
    star_index: Path,
    config: Optional[InferenceConfig] = None,
) -> Tuple[StructureInference, Tuple[JunctionCheck, JunctionCheck]]:
    """Align the emitted inserts and fold the verdicts into Q2/Q3, the
    transform decision and the architecture."""
    config = config or InferenceConfig()
    if result.architecture is None:
        raise ValueError("no architecture to check")
    profile = run_star(
        emitted_inserts(result.architecture, reads, config), star_index, config
    )
    q2 = check_junction(result.q2, "5", profile, config)
    q3 = check_junction(result.q3, "3", profile, config)
    junctions = replace(result.junctions, q2=q2.answer, q3=q3.answer)
    transform = decide_transform(
        result.observation, result.q1, q2.answer, q3.answer, config
    )
    transform = replace(
        transform, flags=transform.flags + tuple(c.flag for c in (q2, q3) if c.flag)
    )
    architecture = build_architecture(result.q1, junctions, transform, config)
    checked = replace(
        result, junctions=junctions, transform=transform, architecture=architecture
    )
    return checked, (q2, q3)
