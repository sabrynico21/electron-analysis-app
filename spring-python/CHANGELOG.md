# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.0] - 2026-10-20

First scientific release. The package is stable and its numerical behaviour
against the R reference is documented in `docs/EQUIVALENCE.md`.

### Added
- Initial port of SPRING (`SPRING.default`), `hugeKmb`, `mclr`.
- Faithful re-implementations of `latentcor` (truncated-type latent
  correlation), `huge`'s MB coordinate descent (2.0.1), `Matrix::nearPD`, and
  `pulsar`'s StARS.
- Shipped `ipol_22` interpolation table extracted from `latentcor` 2.0.2.
- CLI (`python -m spring_microbiome`) with a JSON stdin→stdout contract for the
  Electron host.
- Golden R fixtures and the equivalence test-suite.
- Golden coverage at the **SPRING defaults** (`nlambda=20`, `rep.num=20`,
  `thresh=0.1`) for both quantitative and compositional (`mclr`) input.
- Golden coverage of `Rmethod="original"` (numeric bridge inversion) against R
  through a dedicated latentcor-stage fixture.
- Data-specific λ-path goldens for compositional public data (American Gut and
  GlobalPatterns), exercising the full `latentcor → λ-max → grid` chain.
- Self-contained golden regeneration: the SPRING R source, `QMP.rda` and the
  public-data selector are vendored under `tests/r_reference/reference/`.

### Changed
- `spring-python/` is now tracked by the workspace Git repository; golden
  fixtures are committed so the test-suite never needs R.
- README explicitly declares the divergence from the SPRING paper
  (`latentcor` fork vs. the paper's `mixedCCA`).

[1.0.0]: https://example.invalid/spring-microbiome/releases/tag/v1.0.0
[Unreleased]: https://example.invalid/spring-microbiome/compare/v1.0.0...HEAD
