# Shared

## MS_PeptideEstimator
Takes an input fasta file and chosen proteases, and returns figures displaying the number of cleavage points and the sizes of the smallest possible fragments of proteins.

## BootstrapConvergence
Takes isoform expression feature counts from Kallisto, and finds the isoforms which have the largest difference between deterministic feature counts and bootstrap estimates which are used in downstream software for technical variance estimations. This was used to decide on an appropriate number of bootstraps to perform on data.
