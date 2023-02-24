# smunger


[![pypi](https://img.shields.io/pypi/v/smunger.svg)](https://pypi.org/project/smunger/)
[![python](https://img.shields.io/pypi/pyversions/smunger.svg)](https://pypi.org/project/smunger/)
<!-- [![Build Status](https://github.com/jianhua/smunger/actions/workflows/dev.yml/badge.svg)](https://github.com/jianhua/smunger/actions/workflows/dev.yml) -->
<!-- [![codecov](https://codecov.io/gh/jianhua/smunger/branch/main/graphs/badge.svg)](https://codecov.io/github/jianhua/smunger) -->
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)



munger for GWAS summary statistics


<!-- * Documentation: <https://jianhua.github.io/smunger> -->
<!-- * GitHub: <https://github.com/jianhua/smunger> -->
* PyPI: <https://pypi.org/project/smunger/>
* Free software: MIT


## Features

- [x]  define column properties
    - [x]  required columns: CHR, BP, EA, NEA
    - [x]  optional columns: BETA, SE, P, EAF, MAF
    - [x]  Auxiliary columns: OR, OR_SE, Z
    - [x]  Data types
    - [x]  Data ranges
    - [x]  Allow missing values and default missing values
- [x]  semi-automatically header mapping
    - [x]  read first five rows and display in terminal
    - [x]  guess header map by common column names
    - [x]  manually check if the mapping is correct
    - [x]  input the right column number if it is wrong
    - [x]  check if OR, OR_SE, Z are present if BETA, SE are absent
    - [x]  save the final column map to json for further munging
- [x]  data munging
    - [x]  EA ≠ NEA
    - [x]  if EAF presents, MAF = min(EAF, 1-EAF)
    - [x]  convert OR/ORSE to BETA/SE, if BETA, SE are absent and OR, ORSE are present
    - [x]  remove duplicate SNPs with same chr-bp-sorted(EA,NEA), keep the one with lowest P
    - [x]  output: \t separated, `bgzip` compress, `tabix` index.
    - [x]  optional output: significant SNPs, munge report
    
    |  | CHR | BP | rsID | EA | NEA | EAF | MAF | BETA | SE | P | OR | OR_SE | Z |
    | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
    | type | int | int | str | str | str | float | float | float | float | float | float | float | float |
    | allow null | False | False | True | False | False | False | False | True | False | True | True | False | True |
    | null value |  |  |  |  |  |  |  | 0 |  | 0.999 | 1 |  | 0 |
    | range | [1，23] | (0,inf) |  | only contains ‘ACGT’ | only contains ‘ACGT’ | [0,1] | [0,0.5] | (-inf,inf) | (0, inf) | (0,1) | (0, inf) | (0, inf) | (-inf,inf) |
- [x]  liftover
    - [x]  guess genome build
    - [x]  liftover
- [x]  annotate
    - [x]  annotate rsID

## Credits

This package was created with [Cookiecutter](https://github.com/audreyr/cookiecutter) and the [waynerv/cookiecutter-pypackage](https://github.com/waynerv/cookiecutter-pypackage) project template.
