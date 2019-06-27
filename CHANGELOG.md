# Changelog

## [v1.0.0rc3-1]

### Fixed
- Fixed bug in CLI entrypoint

## [v1.0.0rc3]

### Added
- Reads simulator now uses variability measured from the data
- Switched to Poetry
- Improved error reporting in TxComp

### Fixed
- Fixed multithreading errors
- Fixed error in SampCompDB.plot_position() (#85)
- Fixed errors with 0 pvalues (#87 and #90)
- Fixed error when passing a Whitelist object to SampComp (#91)

## [v1.0.0rc2]

### Added
- Continuous testing with Travis CI
- Automatic deployment of docs to gh-pages

### Fixed
- Fixed "Not enough p-values" error. Issue #68
