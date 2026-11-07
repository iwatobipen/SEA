# Calculation method and differences from the original paper

1. Standardize structures using RDKit Cleanup → FragmentParent → Uncharger → canonical isomeric SMILES. Remove duplicate structures within each set and retain shared structures across sets. Tautomers are not merged.
2. Generate RDKit path fingerprints with 2048 bits, `minPath=1`, and `maxPath=7`. Other settings use the defaults of the locked RDKit version. The raw score is the sum of Tanimoto coefficients for all molecule pairs at or above the threshold; the threshold is not subtracted from each contribution.
3. Sample calibration molecules uniformly from the unique structures in the reference database. Sample without replacement within each set and independently between sets, allowing structures to overlap across sets.
4. Generate repeated samples for each size pair. Use the first half to estimate the mean `mu(n)=a*n` and standard deviation `sigma(n)=b*n**c`, where `n` is the product of the two set sizes. Fit the standard deviation in log space.
5. Calculate Z-scores and fit Gumbel location and scale parameters for thresholds from 0.30 to 0.80 in increments of 0.01. Select the threshold with the lowest KS statistic on the second half of the samples. Exclude size pairs with fewer than 10 molecules in either set or zero score variance from distribution fitting. Use a separate scaffold holdout for final target-retrieval evaluation.
6. Use Gumbel p-values when the held-out KS statistic is ≤ 0.10, both sets have at least 10 molecules, and their sizes are within the calibrated range. Otherwise, use Monte Carlo evaluation at the actual set sizes. The default is 999 trials, with `p=(exceed+1)/(trials+1)`. Reuse background scores for matching size pairs across jobs in the same worker, isolated by dataset ID and storage directory. A zero raw score has exactly `p=1` because all background raw scores are nonnegative; return `method=exact_zero` without sampling in that case. Reject inputs whose set sizes exceed the background pool size.
7. Calculate `E=p*N`, where `N` is the total number of searched targets, regardless of `top_k`. For a two-set comparison, `N=1`. Log p-values and log E-values use natural logarithms. E-values may exceed 1. Rank results by ascending log E-value, breaking ties by descending raw score. Finite Monte Carlo resolution can produce ties.

These are implementation choices within the original paper's statistical framework; they do not guarantee the same background distribution. In particular, single-molecule queries use finite-sample Monte Carlo evaluation instead of assuming a continuous extreme-value distribution. Threshold-selection and Gumbel-fit diagnostics are stored in the dataset's model.

The original paper's coefficients are used only in dedicated tests. The numerical example in the original Supplementary Methods contains rounding differences, so exact agreement with the E-value reported in the paper is not expected.

- Original paper: https://doi.org/10.1038/nbt1284
- Supplementary Methods: https://static-content.springer.com/esm/art%3A10.1038%2Fnbt1284/MediaObjects/41587_2007_BFnbt1284_MOESM64_ESM.pdf
- ChEMBL assay confidence: https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/chembl-data-questions
- RDKit: https://www.rdkit.org/docs/GettingStartedInPython.html

The initial release does not include all-target networks, an MCP server, public-facing account management, or distributed workers. The single-host SQLite queue is intended for small shared lab deployments.
