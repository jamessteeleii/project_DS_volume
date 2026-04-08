# A test of higher and lower fractional volumes of resistance training upon arm and thigh muscle area: A multi-site randomised trial

## Abstract
Recent work has theorised the effects of resistance training volume to be positive and monotonic, albeit with diminishing returns, with regards to hypertrophy. Improvements in muscle size however are typically small, even smaller in trained people due to the linear-logarithmic adaptation to RT over time, and thus between intervention differences in effects are likely to be very small. As such, in contrast to most studies in the field which aim to detect differences between interventions, we sought to conduct a highly powered pre-registered test of the statistical equivalence of two RT interventions in previously trained participants; namely low (9 fractional sets per week) and high (36 fractional sets per week) volumes. A randomised controlled trial across 22 sites was employed with 125 partcipants recruited. Our primary outcome was hypertrophy operationalised as estimated muscle cross sectional area using circumference and skinfold measurements of the upper arm and thigh. At the participant level, 120 participants were randomly assigned to either the low (n = 56) or high (n = 64) volume RT intervention condition. Participants underwent pre-intervention testing and then participated in a 12-week intervention with post-intervention testing following this. Our primary estimand of interest was the condition by time interaction effect from our pre-registered analysis of pooled outcomes reflecting the standardised between condition difference in change in hypertrophy over time. After randomisation 112 participants completed baseline testing and 87 completed post-intervention testing. The estimate for this effect was 0.023 [95%CI: -0.044, 0.091] and the *p*-value for equivalence was *p*=0.032 supporting statistically equivalent effects between conditions. Main effects for time were also small 0.087 [95%CI: 0.053, 0.121] in line with prior predictions from theoretical linear-log growth models. This study is to our knowledge one of the largest to compare the effects of low and high volume RT interventions upon hypertrophy in previously trained participants. We found statistical equivalence between both conditions and both main effects of time, and any interaction effects for condition by time, are likely small. More broadly, this study further corroborates the linear-log theory of adaptation, that the effects of RT in trained persons should be expected to be small, and that current studies in the field of RT are woefully underpowered to be able to detect their effects, let alone test between intervention comparisons.

## Reproducibility
This repository contains the necessary files and code to reproduce the analyses, figures, and the manuscript. 

## Usage
To reproduce the analyses, you will need to have R (https://cran.r-project.org/) and RStudio (https://www.rstudio.com/products/rstudio/download/#download) installed on your computer.

To help with reproducibility, this project uses the `renv` R package (see https://rstudio.github.io/renv/articles/renv.html). With `renv`, the state of this R project can be easily loaded as `renv` keeps track of the required R packages (including version), and (if known) the external source from which packages were retrieved (e.g., CRAN, Github). With `renv`, packages are installed to a project specific library rather than your user or system library. The `renv` package must be installed on your machine before being able to benefit from its features. The package can be installed using the following command:

``` r
install.packages("renv")
```

Once you have `renv` installed, you can get a copy of this repository on your machine by clicking the green Code button then choose Download zip. Save to your machine and extract. After extraction, double click the `project_DS_volume.Rproj` file in the root directory. This will automatically open RStudio. This will ensure all paths work on your system as the working directory will be set to the location of the `.Rproj` file. Upon opening, RStudio will recognize the `renv` files and you will be informed that the project library is out of sync with the lockfile. At shown in the console pane of RStudio, running `renv::restore()` will install the packages recorded in the lockfile. This could take some time depending on your machine and internet connection.

## Targets analysis pipeline

This project also uses a function based analysis pipeline using
[`targets`](https://books.ropensci.org/targets/). Instead of script based pipelines the `targets` package makes use of functions applied to targets specified within the pipeline. The targets can be viewed in the `_targets.R` file, and any user defined functions are available in `R/functions.r`.

You can view the existing targets pipeline by clicking [here](https://jamessteeleii.github.io/project_DS_volume/targets_pipeline.html).

Useful console functions:

- `tar_edit()` opens the make file
- `tar_make()` to run targets
- `tar_visnetwork()` to view pipeline

## Software and packages used

The [`grateful`](https://pakillo.github.io/grateful/index.html) package was used to create citations to all software and packages used in the analysis. The `grateful` report can be viewed by downloading by clicking [here](https://jamessteeleii.github.io/project_DS_volume/grateful-report.html).

## License

Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
  [cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg

