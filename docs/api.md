# API Reference

## Core Functions

### Cleanliness Checking

```python
from getRPF.core import handle_cleanliness_check

handle_cleanliness_check(
    input_file="reads.fastq",
    format="fastq",
    output="report.txt",
    min_quality=20,
    threads=1,
    rpf_checks=True
)
```

### Adapter Detection

```python
from getRPF.core import handle_adapter_detection

handle_adapter_detection(
    input_file="reads.fastq",
    format="fastq",
    output="report.txt",
    adapter="AGATCGGAAGAG",
    min_overlap=10,
    max_mismatches=1,
    threads=1
)
```

## Check Classes

### Length Distribution Check

```python
from getRPF.core.checkers import LengthDistributionCheck

checker = LengthDistributionCheck(
    min_length=26,
    max_length=35,
    warn_threshold=0.5,
    fail_threshold=0.3
)
```

### Base Composition Check

```python
from getRPF.core.checkers import BaseCompositionCheck

checker = BaseCompositionCheck(
    max_base_freq=0.85
)
```

### GC Content Check

```python
from getRPF.core.checkers import GCContentCheck

checker = GCContentCheck(
    min_gc=0.35,
    max_gc=0.65
)
```

## Result Classes

### CleanlinessResults

```python
from getRPF.core.processors.check import CleanlinessResults

results = CleanlinessResults(
    length_distribution={},
    nucleotide_frequencies={},
    quality_scores=None,
    gc_content=None
)
```

### CheckResult

```python
from getRPF.core.checkers import CheckResult, Status

result = CheckResult(
    status=Status.PASS,
    message="Check passed successfully",
    details={}
)
```
