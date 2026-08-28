# Racial Gaps in Managerial Career Advancement

Replication package for the paper by Sam Hoey, Thomas Peeters and Stefan Szymanski.

**Full documentation, including the data dictionary, is in [README.pdf](README.pdf).**
This page is a summary.

## What this is

Code and data to reproduce every table and figure in the paper. The package contains four R
scripts, the pseudonymised data they read, and the output they produce.

The data come from Transfermarkt and from a purpose built picture survey used to code race.
Player, coach and club names, profile URLs and exact dates of birth are replaced by stable
pseudonyms, so the distributed files carry no direct identifiers. Only the year of birth is
retained, because the specifications use year of birth fixed effects.

## Contents

| Path | Holds |
|---|---|
| `HPS Data preparation.R` | Builds the analysis environment from the raw inputs |
| `HPS Simulation preparation.R` | Fits the transition model and runs the simulations |
| `HPS Racial Gaps.R` | Every result except the simulation ones |
| `HPS Simulation results.R` | The simulation results |
| `data/source/` | Raw inputs, read only by the preparation code |
| `data/` | The environments the preparation code builds |
| `output/tables/` | 13 tables, as `.tex` |
| `output/figures/` | 26 numbered figures, as `.pdf` |
| `output/other/` | 5 figures the code produces that the paper does not report |

Output files are named for the table or figure they correspond to in the paper, so
`output/tables/table4_hazard_model.tex` is Table 4 and `output/figures/figure6_survival_2x2.pdf`
is Figure 6.

## Running it

Set the working directory to the package root, then

```r
setwd("/path/to/this/folder")
source("HPS Data preparation.R")
source("HPS Simulation preparation.R")   # once for each uksim setting
source("HPS Racial Gaps.R")
source("HPS Simulation results.R")
```

The simulation environments are included, so all results can be reproduced without repeating
the simulations, which take several hours per sample.

Required packages: `collapse`, `plyr`, `dplyr`, `tidyverse`, `ggplot2`, `lfe`, `xtable`,
`RColorBrewer`, `survival`, `Matching`, `nnet`, `data.table`, `stringr`.

## Licence

The R scripts are released under the MIT License. The data and output are released under
CC BY 4.0. See [LICENSE](LICENSE) for the full terms, including a note on the third party
provenance of the underlying match and career records.
