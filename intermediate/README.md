## Intermediate Result Summaries
This `intermediate` subdirectory contains three subdirectories (`gard`, `iqtree`, and `beast`), each corresponding to one type of analysis and comprising the results summarized from the posterior outputs (which can then used as inputs for making the result figures and tables).
The files provided in this subdirectory intend to facilitate a quick reproduction of the final result figures and tables presented in our manuscript and supplementary material.

The `gard` subdirectory contains a spreadsheet showing the AIC-c score computed by `GARD` under various numbers of breakpoints.
Each row corresponds to a dataset (11 rows total excluding the first header row as `GARD` was unable to run for four of the 15 datasets because their alignments are shorter than what `GARD` requires).
The recombination scheme with three breakpoints was only evaluated for the rabies virus dataset.
The `GARD` analyses of two datasets (`ebov_dud17` and `sars2_lem21`) failed to complete within the 30-day cluster walltime so the AIC-c scores are only available for the recombination free scheme.
The AIC-c score is recorded as `NA` for the recombination schemes that were not evaluated.

Each dataset-specific subdirectory in `iqtree` contains the following four outputs summarized from the IQ-TREE analyses and the subsequent RTT regression test :
| file name | variables |
|----------|----------|
| mlscore.tsv    | the (log) maximum likelihood (ML) value of each replicated IQ-TREE run|
| ml.tree   | the ML tree inferred across all replicate runs |
| ml_rootbyrtt.tree   | the ML tree rooted by RTT regression |
| rttvals.rds   | results of the RTT regression test |

Each dataset-specific subdirectory in `beast` contains the MCMC diagnostics computed with the posterior samples in the corresponding `analyses` subdirectory.
The `ess_*.tsv` and `psrf_*.tsv` spreadsheet files present the effective sample size (ESS) and potential scale reduction factor (PSRF) of each major variable, respectively; similarly, the `treess_*.tsv` and `treepsrf_*.tsv` store the tree ESSs and PSRFs.
The `pairwiseR*2d.rds` contains the MDS coordinates of each posterior tree sample in the reduced two-dimensional space.