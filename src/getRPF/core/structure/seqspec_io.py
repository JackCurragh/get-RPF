"""Architecture <-> seqspec YAML, and catalogue templates as hypotheses.

Seqspec is the file format; the Architecture is the model (spec §9). Files
written here carry a ``getrpf`` key per region, so they round-trip exactly.
Plain seqspec files, and the catalogue in ``getRPF/architectures``, are read
by region type. Templates are compared with an inferred architecture and
marked supported, partially supported or rejected (spec §5.5); they never
decide a boundary.
"""

from __future__ import annotations

import importlib.resources
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Sequence, Tuple, cast

import yaml

from .model import Architecture, Block, Status, TransformDecision

_WRITE_TYPE = {
    "random": "umi",
    "umi": "umi",
    "barcode": "barcode",
    "fixed": "linker",
    "nta": "linker",
    "insert": "cdna",
    "adapter": "adapter",
}
_SEQUENCE_TYPE = {
    "random": "random",
    "umi": "random",
    "barcode": "random",
    "fixed": "fixed",
    "nta": "random",
    "insert": "joined",
    "adapter": "fixed",
    "tail": "fixed",
}


class _Loader(yaml.SafeLoader):
    """Safe loader that reads unknown tags (``!Region``, python tags) as data."""


def _untagged(loader: yaml.SafeLoader, suffix: str, node: yaml.Node) -> Any:
    if isinstance(node, yaml.MappingNode):
        return loader.construct_mapping(node)
    if isinstance(node, yaml.SequenceNode):
        return loader.construct_sequence(node)
    return loader.construct_scalar(cast(yaml.ScalarNode, node))


_Loader.add_multi_constructor("!", _untagged)
_Loader.add_multi_constructor("tag:yaml.org,2002:python/", _untagged)


def to_seqspec(
    architecture: Architecture,
    transform: Optional[TransformDecision] = None,
    assay_id: str = "getrpf_inferred",
) -> Dict[str, Any]:
    regions = []
    for index, block in enumerate(architecture.blocks):
        if block.type == "tail":
            region_type = f"poly_{block.sequence}"
        elif block.type == "random" and not block.keep_as_umi:
            region_type = "linker"
        else:
            region_type = _WRITE_TYPE[block.type]
        region: Dict[str, Any] = {
            "region_id": f"{block.type}_{index}",
            "region_type": region_type,
            "sequence_type": _SEQUENCE_TYPE[block.type],
            "min_len": block.length[0],
            "max_len": block.length[1],
        }
        if block.sequence:
            region["sequence"] = block.sequence
        region["getrpf"] = {
            "type": block.type,
            "frame": block.frame,
            "keep_as_umi": block.keep_as_umi,
            "remove": block.remove,
            "status": block.status.value,
            "note": block.note,
        }
        regions.append(region)
    meta: Dict[str, Any] = {
        "source": architecture.source,
        "status": architecture.status.value,
        "fragment_policy": architecture.fragment_policy,
        "describe": architecture.describe(),
    }
    if transform is not None:
        meta["transform"] = {
            "emit": transform.emit,
            "bound": transform.bound.value,
            "reasons": list(transform.reasons),
            "conventions": list(transform.conventions),
            "flags": list(transform.flags),
        }
    return {
        "seqspec_version": "0.3.0",
        "assay_id": assay_id,
        "name": "get-RPF inferred read structure",
        "modalities": ["rna"],
        "sequence_spec": regions,
        "getrpf": meta,
    }


def write_seqspec(
    path: Path,
    architecture: Architecture,
    transform: Optional[TransformDecision] = None,
) -> None:
    path.write_text(
        yaml.safe_dump(to_seqspec(architecture, transform), sort_keys=False)
    )


def read_seqspec(path: Path) -> Tuple[Architecture, Optional[Dict[str, Any]]]:
    """The architecture, and the transform decision if infer-structure wrote one."""
    spec = load_yaml(path.read_text())
    meta = spec.get("getrpf") or {}
    return from_seqspec(spec), meta.get("transform")


def load_yaml(content: str) -> Dict[str, Any]:
    data = yaml.load(content, Loader=_Loader)  # noqa: S506 - safe loader subclass
    if not isinstance(data, dict):
        raise ValueError("seqspec is not a mapping")
    return data


def from_seqspec(spec: Dict[str, Any]) -> Architecture:
    regions = list(
        _flatten(spec.get("sequence_spec") or spec.get("library_spec") or [])
    )
    at = next(
        (
            i
            for i, region in enumerate(regions)
            if str(region.get("region_type", "")).lower() in ("cdna", "rpf")
        ),
        None,
    )
    if at is None:
        raise ValueError("seqspec has no cdna/rpf region")
    blocks = tuple(
        _block(region, "read_start" if i <= at else "anchor", i == at)
        for i, region in enumerate(regions)
    )
    meta = spec.get("getrpf") or {}
    low, high = blocks[at].length
    return Architecture(
        blocks,
        cast(Any, meta.get("source", "template")),
        Status(meta.get("status", "resolved")),
        meta.get("fragment_policy") or f"template_{low}_{high}",
    )


def _flatten(regions: Sequence[Any]) -> Iterator[Dict[str, Any]]:
    for region in regions:
        if not isinstance(region, dict):
            continue
        if region.get("regions"):
            yield from _flatten(region["regions"])
        else:
            yield region


def _block(region: Dict[str, Any], frame: str, is_insert: bool) -> Block:
    low = int(region.get("min_len") or 0)
    high = int(region.get("max_len") or low)
    sequence = region.get("sequence") or next(iter(region.get("sequences") or []), None)
    meta = region.get("getrpf")
    if meta:
        return Block(
            cast(Any, meta["type"]),
            cast(Any, meta["frame"]),
            (low, high),
            sequence,
            bool(meta["keep_as_umi"]),
            bool(meta["remove"]),
            Status(meta["status"]),
            meta.get("note", ""),
        )
    region_type = str(region.get("region_type", "")).lower()
    sequence_type = str(region.get("sequence_type", "")).lower()
    resolved = Status.RESOLVED
    edge = cast(Any, frame)
    if is_insert:
        return Block("insert", "read_start", (low, high), None, False, False, resolved)
    if region_type == "umi":
        return Block("random", edge, (low, high), None, True, True, resolved)
    if region_type == "barcode":
        return Block("barcode", edge, (low, high), sequence, False, True, resolved)
    if region_type == "adapter":
        return Block("adapter", "anchor", (low, high), sequence, False, True, resolved)
    if region_type.startswith("poly_"):
        base = region_type[len("poly_") :].upper()
        return Block("tail", "anchor", (0, high), base, False, True, resolved)
    kind = "fixed" if sequence_type == "fixed" else "random"
    return Block(cast(Any, kind), edge, (low, high), sequence, False, True, resolved)


# --- Templates as hypotheses (spec §5.5) ---


@dataclass(frozen=True)
class TemplateComparison:
    template: str
    verdict: str
    """supported | partially_supported | rejected"""
    reason: str


def load_catalogue() -> List[Tuple[str, Architecture]]:
    """Every parseable template in the package catalogue."""
    templates: List[Tuple[str, Architecture]] = []
    directory = importlib.resources.files("getRPF").joinpath("architectures")
    for entry in sorted(directory.iterdir(), key=lambda item: item.name):
        if not entry.name.endswith((".yaml", ".yml")):
            continue
        try:
            spec = load_yaml(entry.read_text())
            templates.append(
                (str(spec.get("assay_id") or entry.name), from_seqspec(spec))
            )
        except (yaml.YAMLError, ValueError, KeyError, TypeError):
            continue
    return templates


def compare_catalogue(
    architecture: Architecture,
    catalogue: Optional[Sequence[Tuple[str, Architecture]]] = None,
) -> List[TemplateComparison]:
    comparisons = [
        compare(architecture, name, template)
        for name, template in (catalogue if catalogue is not None else load_catalogue())
    ]
    rank = {"supported": 0, "partially_supported": 1, "rejected": 2}
    return sorted(comparisons, key=lambda c: (rank[c.verdict], c.template))


def compare(
    inferred: Architecture, name: str, template: Architecture
) -> TemplateComparison:
    """Compare structure features. Unspecified template features are not held
    against it; a poly(A)-like template 'adapter' is read as a tail."""
    a, b = _features(inferred), _features(template)
    matches: List[str] = []
    mismatches: List[str] = []
    for key, label in (("five", "5' technical"), ("three", "3' technical")):
        if a[key] == b[key]:
            matches.append(f"{label} {a[key]} nt")
        else:
            mismatches.append(f"{label} {a[key]} nt vs {b[key]} nt")
    if a["tail"] == b["tail"]:
        if a["tail"]:
            matches.append(f"poly({a['tail']}) tail")
    else:
        mismatches.append(f"tail {a['tail'] or 'none'} vs {b['tail'] or 'none'}")
    if a["adapter"] and b["adapter"]:
        if _shares_kmer(a["adapter"], b["adapter"]):
            matches.append("adapter family")
        else:
            mismatches.append("adapter differs")
    if not mismatches:
        verdict = "supported"
    elif len(matches) >= 2:
        verdict = "partially_supported"
    else:
        verdict = "rejected"
    reason = "; ".join(
        part
        for part in (
            ("matches: " + ", ".join(matches)) if matches else "",
            ("differs: " + ", ".join(mismatches)) if mismatches else "",
        )
        if part
    )
    return TemplateComparison(name, verdict, reason)


def _features(architecture: Architecture) -> Dict[str, Any]:
    blocks = architecture.blocks
    at = next(i for i, block in enumerate(blocks) if block.type == "insert")
    adapter = next(
        (b.sequence for b in blocks if b.type == "adapter" and b.sequence), None
    )
    tail = next((b.sequence for b in blocks if b.type == "tail"), None)
    if adapter and adapter.count("A") >= 0.7 * len(adapter):
        adapter, tail = None, "A"
    return {
        "five": sum(b.length[1] for b in blocks[:at] if b.remove),
        "three": sum(
            b.length[1]
            for b in blocks[at + 1 :]
            if b.remove and b.type not in ("adapter", "tail")
        ),
        "adapter": adapter,
        "tail": tail,
    }


def _shares_kmer(a: str, b: str, k: int = 10) -> bool:
    kmers = {a[i : i + k] for i in range(len(a) - k + 1)}
    return any(b[i : i + k] in kmers for i in range(len(b) - k + 1))
