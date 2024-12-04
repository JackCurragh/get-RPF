# get-RPF
# getRPF

Get Ribosome Protected Fragment features - a modular toolkit for analyzing Ribo-seq read composition

## Project Structure

```
getRPF/
├── src/
│   └── getRPF/
│       ├── __init__.py          # Package initialization, version info
│       ├── core/                # Core functionality and base classes
│       │   ├── __init__.py
│       │   └── reader.py        # Abstract base classes for file reading
│       ├── io/                  # Input/Output handling
│       │   ├── __init__.py
│       │   ├── fasta.py         # FASTA format handlers
│       │   ├── fastq.py         # FASTQ format handlers
│       │   └── compression.py   # Compression handling (gzip)
│       ├── quality/             # Quality assessment modules
│       │   ├── __init__.py
│       │   └── metrics.py       # Quality metrics calculations
│       ├── detectors/           # Feature detection modules
│       │   ├── __init__.py
│       │   ├── adapter.py       # Adapter detection
│       │   ├── umi.py          # UMI detection
│       │   └── barcode.py      # Barcode detection
│       └── utils/              # Utility functions
│           ├── __init__.py
│           └── helpers.py      # Common helper functions
├── tests/                     # Test directory
│   ├── __init__.py
│   ├── test_io/              # IO tests
│   ├── test_quality/         # Quality module tests
│   └── test_detectors/       # Detector tests
├── docs/                     # Documentation
│   ├── conf.py              # Sphinx configuration
│   └── index.rst            # Documentation root
├── examples/                 # Example usage scripts
├── README.md                # This file
├── pyproject.toml           # Project metadata and build configuration
├── setup.cfg               # Package configuration
└── .gitignore             # Git ignore rules
```

## Development Environment

### Python Version
- Python 3.9 or higher required
- Type hints and f-strings heavily used
- Async IO support for efficient file handling

### Core Dependencies
- `biopython`: Biological sequence parsing
- `numpy`: Numerical operations
- `pandas`: Data manipulation
- `pydantic`: Data validation and settings management

### Development Tools

#### Code Quality
- **Type Checking**: `mypy`
  - Strict type checking enabled
  - Custom type stubs for third-party packages
- **Formatting**: `black`
  - Line length: 88 characters
  - Skip string normalization
- **Linting**: `ruff`
  - Combines multiple linters
  - Autofix capability
  - Strict settings enabled
- **Import Sorting**: `isort`
  - Compatible with Black
  - Separate sections for stdlib, third-party, and local imports

#### Testing
- **Framework**: `pytest`
  - Fixture-based testing
  - Parametrized tests for multiple scenarios
- **Coverage**: `pytest-cov`
  - Minimum coverage requirement: 90%
- **Property Testing**: `hypothesis`
  - For complex input scenarios
  - Custom strategies for biological sequences

#### Documentation
- **Generator**: `sphinx`
  - Google-style docstrings
  - Auto-generated API documentation
- **Format**: `sphinx-rtd-theme`
  - ReadTheDocs theme
  - Mobile-friendly layout

### Development Setup

1. Clone the repository:
```bash
git clone https://github.com/yourusername/getRPF.git
cd getRPF
```

2. Create a virtual environment:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: .\venv\Scripts\activate
```

3. Install development dependencies:
```bash
pip install -e ".[dev]"
```

4. Install pre-commit hooks:
```bash
pre-commit install
```

### Running Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=getRPF

# Run specific test files
pytest tests/test_io/
```

### Building Documentation

```bash
# Generate HTML documentation
cd docs
make html
```

### Code Quality Checks

```bash
# Run type checking
mypy src/getRPF

# Run linter
ruff check src/

# Format code
black src/
isort src/
```

## Contributing

Please read our [Contributing Guidelines](CONTRIBUTING.md) before submitting pull requests.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
