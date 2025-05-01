Source code for reproducing model runs, tables and figures of the manuscript: 
"Formation of spatial vegetation patterns in heterogeneous environments" by
Karl Kästner et al. submitted to Plos One.

Requirement:
- GNU/Linux (recommended operating system)
- Matlab R2023b (older versions probably work)
- dependencies from https://github.com/karlkastner
  (These are fetched automatically when running the batch script)

The script "pattern_formation_batch" reproduces the whole workflow.
It first fetches dependencies from GitHub, and subsequently executes scripts 
for fetching and analyzing the patterns, followed by plotting the results.

Note that the configuration for reproducing the data and figures used in the
manuscript require in 200 long simulations and several shorter simulations.
The runtime of a long simulation is about 1 day and of all simulations combined
more than 6 month on a single cpu-core. For testing, the runtime can be
reduced by reducing the number of cases of degrees of heterogeneity cva, 
the domain size L or simulated time T by adapting the respective parameters in
the configuration file rk_2d_heterogeneity_setup.

The GPL v. 3 extends over all files in the repository, i.e. computer scripts,
documentation and geospatial data files, except for the sample satellite images
in the input folder, which are property of Google / Maxar Technologies (2023).  

