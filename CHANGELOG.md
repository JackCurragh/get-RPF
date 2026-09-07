# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Releases before 0.4.0 predate this file. Five tags (v0.2.1 through v0.3.1) were
cut, none of which published to PyPI, and none of which had a test gate in front
of them.

## [Unreleased]

### Added
- `ci.yml`, the single definition of the project gate: ruff + black, mypy,
  tests on Python 3.10 and 3.12, and a build job that asserts both the
  normalised artefact filenames and the presence of the architecture YAMLs in
  the wheel. Until now the test suite ran nowhere in CI.
- `release.yml`, which calls that gate on a `v*` tag and then publishes to
  PyPI via Trusted Publishing, pushes the release container tags, and creates
  a GitHub Release from the artefacts the gate built.
- `.bumpversion.cfg` covering `pyproject.toml`, `src/getRPF/__init__.py` and
  `CITATION.cff`; the version was previously duplicated with nothing syncing
  the copies.
- `CITATION.cff`, `CHANGELOG.md`, `Makefile` and `mypy.ini`.

### Changed
- Build backend is setuptools rather than hatchling, and the distribution name
  is `getrpf` rather than `getRPF`. The import package, the console script and
  the container entrypoint are unchanged.
- Container tags now have one owner each: `:latest`, `:vX.Y.Z` and `:X.Y` come
  from a release, `:main` and `:sha-<short>` from the tip of `main`. The
  workflow that wrote both, with no test gate, is deleted.
- `isort` and flake8-era configuration removed; ruff covers import order.

### Fixed
- Soft-clip analysis raised `TypeError` on an aligned record with no SEQ
  (`query_sequence` is `Optional`); such records are now skipped with the
  other unusable ones.
- `seqspec_loader` could derive an RPF length range of `(None, None)` when a
  seqspec region omitted `min_len`/`max_len`. Unbounded regions are skipped
  and the documented `(20, 40)` default applies.
- `_find_best_match` call sites narrowed on the match position alone and then
  indexed the partial-match histogram with a possibly-`None` match length.
- `handle_cleanliness_check` and `handle_adapter_check` forwarded `None`
  through to processor parameters typed as `int`; they now fall back to the
  defaults those processors declare.
- `MIN_MODE_FREQUENCY` in the soft-clip analyser was defined and then ignored
  in favour of an inline literal.
