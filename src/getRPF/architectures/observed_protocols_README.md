# Observed Protocol Architectures

These seqspec YAML files encode adapter and UMI structures observed across public
Ribo-seq runs. They are bundled as generic protocol-shape candidates for getRPF
matching and adapter-evidence scoring.

They do not encode accession-level expectations, external database validation,
or dynamic RPF boundaries. Each observed protocol uses a broad RPF length range
so that final output quality is decided from extracted read evidence rather than
from the source metadata that suggested the adapter structure.

Adapter strings are normalized to the DNA alphabet (`U` -> `T`) for matching
sequencing reads.
