# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/althonos/pyncbitk/compare/v0.1.0-alpha.4...HEAD

## [v0.1.0-alpha.5] - 2025-11-04

### Added
- Hash implementation for all `Serial` objects based on ASN serialization.
- `__rich_repr__` implementations for several classes.
- Type hints to `pyncbitk.object.seqset`.
- `descriptions` and `length` properties to `BioSeqHandle` proxy class.
- *WIP* extension `pyncbitk.objects.taxon` providing a client to the NCBI Taxonomy service.

### Changed
- Pin required `pyncbitk-runtime` to `29.6` and later.
- Make `pyncbitk.algo` a package and move BLAST API to `pyncbitk.algo.blast`.

### Fixed
- Typo for MacOS `cibuildwheel` configuration in `pyproject.toml`
- Type name changes in `pyncbitk.objects.seq`.


## [v0.1.0-alpha.4] - 2025-10-30

### Added
- Flag to skip parsing FASTA record IDs in `FASTAReader`.

### Changed
- Reorganized wheel distribution to redistribute `ncbi-cxx-toolkit-public` runtime libraries in a separate `pyncbitk-runtime` wheel.
- Removed unused Conan configuration file from main `pyncbitk` package.
- Update `readthedocs.yaml` configuration to build the package from source on new commits.

### Fixed
- Compilation issues on MacOS related to `$RPATH` management for Python extensions.
- Uninitialized variable in `AlignMap.__getitem__` causing crashes in debug mode only.


## [v0.1.0-alpha.3] - 2025-10-28
[v0.1.0-alpha.3]: https://github.com/althonos/pyncbitk/compare/v0.1.0-alpha.2...v0.1.0-alpha.3

### Added
- Setters for `percent_identity` and `culling_limit` to `Blast`.
- Getters for gaps in `SeqAlign` and individual `AlignRow`.
- `GeneralId.db_tag` property to access the `DBTag` of a `GeneralId`.

### Changed
- Update `cibuildwheel` configuration to strip wheels before uploading to PyPI.


## [v0.1.0-alpha.2] - 2025-10-26
[v0.1.0-alpha.2]: https://github.com/althonos/pyncbitk/compare/v0.1.0-alpha.1...v0.1.0-alpha.2

### Added
- Buffer protocol support to `NcbiEAaData`.
- Encoding constructor to all `SeqData` subclasses.
- Options to configure the `BlastP` word size and word threshold.
- New general class `DBTag` storing a database cross-reference.
- New `pyncbitk.objects.seqdesc` module with some sequence description classes.
- `descriptions` keyword argument to `BioSeq` constructor allowing to set the sequence descriptions.

### Changed
- Make `Blast.run` copy the current options to avoid race conditions.
- Bump NCBI C++ Toolkit to `v29.6.0`.
- Improve NCBIptb build system integration to support editable installs.
- Update Docker image with build artifacts.
- Build wheels in Limited API for Python 3.11 and later.

### Fixed
- Issue in `Blast.run` with arbitrary query iterable not setting the query vector.
- Issue in `Blast.run` on empty subject vector.
- Compilation issues with Conan.
- `DatabaseReader.keys()` iteration not incrementing the underlying C++ iterator.
- Recent versions of `scikit-build-core-conan` not supporting Python 3.8.


## [v0.1.0-alpha.1] - 2024-12-12
[v0.1.0-alpha.1]: https://github.com/althonos/pyncbitk/compare/6fc81aa8...v0.1.0-alpha.1

Initial release.
