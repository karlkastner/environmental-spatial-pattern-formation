Source code for reproducing analysis, tables and figures of the accompanying
manuscript.

Requirement:
- GNU/Linux (recommended operating system)
- Matlab R2023b (older versions probably work)
- dependencies from https://github.com/karlkastner
  (These are fetched automatically when running the batch script)

The script "pattern_formation_batch" reproduces the whole workflow.
It first fetches dependencies from GitHub, and subsequently executes scripts 
for fetching and analyzing the patterns, followed by plotting the results.

Note that the configuration for reproducing the data and figures used in the
manuscript require in total 102 simulations with a runtime of more than 1 day
per simulation on a single cpu-core, in total 3-4 month.
We therefore provide a quick configuration with 12 simulations requiring
each 10 mins, in total 2 hours. The quick configuration is the default.

The GPL v. 3 extends over all files in the repository, i.e. computer scripts,
documentation and geospatial data files, except for the sample satellite images
in the input folder, which are property of Google / Maxar Technologies (2023).  

