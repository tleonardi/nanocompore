# Changelog

## v2.0.0

This is major release that represents an almost complete redesign and rewrite of the tool. The more noteworthy changes include:

- Completely new user interface.
- Support for the RNA004 sequencing chemistry.
- Major improvements in performance.
- GPU support.
- Support for multiple signal-to-reference aligners (resquigglers): Nanopolish/f5c eventalign, Uncalled4, and Remora.
- Reduced disk space usage.
- Automatic test selection strategy.
- We add the ability to use the dwell time of a nearby position as an additional dimension in the GMM test. This can be used to capture changes in the molecule translocation speed produced by interactions between a modification and the motor protein.

## v1.0.4

### Fixed

- Fixed logging levels and verbosity of log messages

## v1.0.3

### Fixed

- Fixed bug Eventalign_collapse CLI options

## v1.0.2

### Added

- Exposed option to enable/disable anova test
- Running nanocompore via CLI now disables anova by default
- SampComp now downsamples to 5000 by default

### Fixed

- Subsampling in whitelist is deterministic (fix for #103)
- Reworked multiprocessing framework for SampComp (fix for zombie threads, not tested with large datasets)

## v1.0.1

### Fixed

- Fixed #120, #122, #138

### Added

- Improved logging

## v1.0.0rc3-1

### Fixed

- Fixed bug in CLI entrypoint

## v1.0.0rc3

### Added

- Reads simulator now uses variability measured from the data
- Switched to Poetry
- Improved error reporting in TxComp

### Fixed

- Fixed multithreading errors
- Fixed error in SampCompDB.plot_position() (#85)
- Fixed errors with 0 pvalues (#87 and #90)
- Fixed error when passing a Whitelist object to SampComp (#91)

## v1.0.0rc2

### Added

- Continuous testing with Travis CI
- Automatic deployment of docs to gh-pages

### Fixed

- Fixed "Not enough p-values" error. Issue #68
