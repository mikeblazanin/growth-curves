# Paper information

## Theoretical validation of growth curves for quantifying phage-bacteria interactions

Michael Blazanin, Emma Vasen, Cèlia Vilaró Jolis, William An, and Paul Turner

**Abstract** Bacteria-infecting viruses, bacteriophages, are the most abundant biological entities on the planet, frequently serving as model systems in basic research and increasingly relevant for medical applications such as phage therapy. A common need is to quantify the infectivity of a phage to a given bacterial host (or the resistance of a host to a phage). However, current methods to quantify infectivity suffer from low-throughput or low-precision. One method that has the potential for high-throughput and high-precision quantification of phage-bacteria interactions is growth curves, where bacterial density is measured over time in the presence and absence of phages. Recent work has proposed several approaches to quantify these curves into a metric of phage infectivity. However, little is known about how these metrics relate to one another or to underlying phage and bacterial traits. To address this gap, we apply ecological modeling of phage and bacterial populations to simulate growth curves across a wide range of trait values. Our findings show that many growth curve metrics provide parallel measures of phage infectivity. Informative metrics include the peak and decline portions of bacterial growth curves, are driven by the interactions between underlying phage and bacterial traits, and correlate with conventional measures of phage fitness. Moreover, we show how intrapopulation trait variation can alter growth curve dynamics. Finally, we test the sensitivity of growth curve metrics to inoculum densities, and assess techniques to compare growth curves across different bacterial hosts. In all, our findings support the use of growth curves for precise high-throughput quantification of phage-bacteria interactions across the microbial sciences.

bioRxiv: 

# Data and Analysis

All of the data and analysis scripts are retained in the github repository https://github.com/mikeblazanin/growth-curves. 

`numerical_analysis5/numerical_analysis5.R` runs all simulations and generates all figures with simulated data.

`Empirical/Empirical_analysis.R` runs all analyses and generates all figures pertaining to experimental data.

# Setup

To re-run the code in this repo, clone the repo then open the associated R Project file. The `renv` package should automatically bootstrap itself, downloading and installing the appropriate version of `renv` into the project library. Once this has completed, run

```
renv::restore()
```

This will install all package dependencies in the version used for our analyses. Then, you can run `numerical_analysis5/numerical_analysis5.R` or `Empirical/Empirical_analysis.R` to generate figures.