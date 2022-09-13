# FrechetLD: Overview

Source code for manuscript:

Frank, S. A. 2022. The number of neutral mutants in an expanding Luria-Delbrück population is approximately Frechet.

by Steven A. Frank, https://stevefrank.org

[License CC-By 4.0](https://creativecommons.org/licenses/by/4.0/)

---

Source code does three things.

1. Takes empirical CDF for a Luria-Delbrück process and fits a Frechet distribution by optimizing the Frechet parameters.
2. Takes an empirical CDF for a Luria-Delbrück process and fits a Frechet distribution by using the parameters given in the manuscript listed above
3. For each case, plot the empirical CDF vs the Frechet fit.

The comments within the file src/Luria.jl describe how to run the code.

The input data for calculating the empirical Luria-Delbrück CDF is at Zenodo at https://doi.org/10.5281/zenodo.7075655. Those data are required to run the code.