.PHONY: help clean clean-build clean-pyc clean-test lint format typecheck test test-cov build dist preflight bump-patch bump-minor bump-major

# Every target runs through one interpreter. Override for a venv or a specific
# version:  make test PYTHON=.pixi/envs/default/bin/python
#
# Tools are invoked as `$(PYTHON) -m <tool>`, never as bare executables: a
# `pip install --user` puts them in a scripts dir that is usually not on PATH.
# The module form is equivalent and PATH-independent (REPO_CONTRACT.md 6).
PYTHON ?= python3

help:
	@echo "PYTHON=$(PYTHON)"
	@echo ""
	@echo "clean        remove build, test and Python artifacts"
	@echo "lint         run ruff + black --check (the CI lint gate)"
	@echo "format       apply black + ruff --fix"
	@echo "typecheck    run mypy (fatal, as in CI)"
	@echo "test         run the test suite"
	@echo "test-cov     run tests with coverage"
	@echo "build/dist   build sdist + wheel"
	@echo "bump-patch   bump patch version, commit and tag (vX.Y.Z)"
	@echo "bump-minor   bump minor version, commit and tag"
	@echo "bump-major   bump major version, commit and tag"
	@echo "preflight    build + twine check + wheel-name + package-data checks"

clean: clean-build clean-pyc clean-test

clean-build:
	rm -fr build/ dist/ .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +

clean-pyc:
	find . -name '*.pyc' -delete
	find . -name '__pycache__' -exec rm -fr {} +

clean-test:
	rm -fr .pytest_cache .mypy_cache .ruff_cache coverage.xml htmlcov/ .coverage

lint:
	$(PYTHON) -m ruff check src tests
	$(PYTHON) -m black --check src tests

format:
	$(PYTHON) -m black src tests
	$(PYTHON) -m ruff check --fix src tests

typecheck:
	$(PYTHON) -m mypy --config-file mypy.ini src/getRPF

test:
	$(PYTHON) -m pytest -q

test-cov:
	$(PYTHON) -m pytest --cov=getRPF --cov-report=term-missing --cov-report=xml

build dist: clean
	$(PYTHON) -m build
	ls -l dist

bump-patch:
	$(PYTHON) -m bumpversion patch

bump-minor:
	$(PYTHON) -m bumpversion minor

bump-major:
	$(PYTHON) -m bumpversion major

# Mirrors ci.yml's `build` job exactly. There is deliberately no TestPyPI
# target: it uses a stored token and exercises neither OIDC nor the
# trusted-publisher binding, so it cannot rehearse the real path
# (REPO_CONTRACT.md 5.4). This is the rehearsal.
preflight: dist
	$(PYTHON) -m twine check dist/*
	@for f in dist/*.whl dist/*.tar.gz; do \
		case "$$(basename $$f)" in \
			getrpf-*) ;; \
			*) echo "ERROR: $$(basename $$f) is not normalised; PyPI needs the 'getrpf-' prefix"; exit 1 ;; \
		esac; \
	done
	@echo "dist filenames normalised:"; ls -1 dist
	@$(PYTHON) -c "import glob, sys, zipfile; \
names = [n for n in zipfile.ZipFile(glob.glob('dist/*.whl')[0]).namelist() \
         if n.startswith('getRPF/architectures/') and n.endswith('.yaml')]; \
print(f'{len(names)} architecture YAMLs in the wheel'); \
sys.exit(0 if names else 'ERROR: no architecture YAMLs in the wheel')"
