# Changelog

## Version 0.3.1

- Make pandas optional. List and search functions now return a `BiocFrame` object.
- Since scipy is only used during upload, the package loads it dynamically and makes it optional.

## Version 0.3.0

- chore: Remove Python 3.8 (EOL).
- precommit: Replace docformatter with ruff's formatter.

## Version 0.2.0

- Changes to support NumPy's 2.0 release.

## Version 0.1.0 - 0.1.3

First release of the package to access, download and save
datasets.
