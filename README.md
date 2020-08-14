# MANA
Modelling Association of Non-Additive effects

## Overview
MANA aims to provide functions to run non-additive GWAS models.

## System Requirements
### Julia Requirements:
  Julia 1.1.0
### Julia dependencies
```
BEDFiles, JWAS
Distributed, DelimitedFiles, Statistics, Linear Algebra, Random, Distributions
DataFrames, JuliaDB, IndexedTables, OnlineStats
```

## Installation Guide
Add package:
```pkg> add https://github.com/egmreynolds/MANA.git```

May need to install BEDFiles
```pkg> add https://github.com/dmbates/BEDFiles.jl.git```

Use package:
```using MANA```

## Example
Minimal example to run non-additive GWAS

``` julia example.jl ```

## License
This project is covered under the MIT License

## Future Updates
Fix issue where initial phenotype LOSO segment does not exist

Add examples for pedigree models

