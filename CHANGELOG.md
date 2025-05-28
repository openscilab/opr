# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
- `delta_h` property 
- `delta_s` property
- Thermodynamic Constants (ΔH, ΔS) calculation (Nearest neighbor)
- Nearest neighbor melting temperature calculation
- `to_protein` method
### Changed
- Test system modified
- `Python 3.6` support dropped
- `README.md` updated
## [0.4] - 2025-03-18
### Added
- Salt-adjusted melting temperature calculation
- `__contains__` overload
- `_computed` internal cache flag attribute
- `is_computed` method
- `to_rna` method
- `E260` property
### Changed
- Cache structure in `Primer` class updated
- Error tests updated
- Cache Test system modified
- `Primer` class `__eq__` method bug fixed
- `README.md` updated
## [0.3] - 2025-02-06
### Added
- `__iter__` overload
- `double_runs` property
- `repeats` method
- `name` property
### Changed
- `single_runs` property updated
- Test system modified
### Removed
- `single_run_length` function
## [0.2] - 2025-01-09
### Added
- `__eq__` overload
- `__str__` overload
- `gc_clamp` property
- `single_runs` property
### Changed
- Test system modified
- `SECURITY.md` updated
- `CONTRIBUTING.md` updated
- `README.md` updated
### Removed
- `property` deleter & setter
## [0.1] - 2024-11-27
### Added
- `MeltingTemperature` enum
- Basic melting temperature calculation
- `addition` and `multipication` operators overload
- `len` magic overload
- `Primer` class
- Molecular weight calculation
- GC content calculation
- Sequence validation
- Unit tests
- `complement` method
- `reverse` method

[Unreleased]: https://github.com/openscilab/opr/compare/v0.4...dev
[0.4]: https://github.com/openscilab/opr/compare/v0.3...v0.4
[0.3]: https://github.com/openscilab/opr/compare/v0.2...v0.3
[0.2]: https://github.com/openscilab/opr/compare/v0.1...v0.2
[0.1]: https://github.com/openscilab/opr/compare/0baa8dd...v0.1