# ChEMBL 37 scoring-parameter validation

Executed on November 8, 2026. The source reference is `e898df53d9bb5240f087`: 1,517 targets and 441,860 unique molecules. The adopted reference ID is `9d52be7ff211a45d2709`. Existing references and completed results remain reproducible.

## Decision and scope

Compare the existing RDKit-path/0.56 configuration with a fixed grid of two representations and 15 thresholds. Select on validation scaffolds, freeze one contender, and assess it on separate test scaffolds. Adopt only if paired 95% bootstrap intervals show an MRR gain and exclude a Recall@10 loss larger than 0.02, at both 32,767 and 199,999 null trials. The production-budget repeat uses the same biological holdout; it is a numerical-resolution check, not an independent second biological replicate.

The selection winner was Morgan radius 2, 2,048 bits, threshold 0.45, without chirality or count simulation. This is a supported improvement within the tested SEA grid, not a universal optimum across representations, assay filters, query types, or future ChEMBL releases.

## What was already ChEMBL-specific

The existing reference did not copy the MDDR coefficients: its threshold was 0.56, mean coefficient 0.0170761465, standard-deviation coefficient 0.0537972549 and exponent 0.7500728844. These were fitted from ChEMBL with 128 samples per size pair. Threshold selection minimized a pooled held-out Gumbel KS statistic over 0.30–0.80; the selected pooled statistic was 0.04166. The selection samples also determined eligibility, and this did not establish per-size or extreme-tail validity.

The 2007 paper used MDDR activity classes and Daylight fingerprints. Dataset and representation differences require new background estimates; a matching numerical Tanimoto threshold does not transfer its statistical meaning.

## Independent null audit

A different seed generated 8,192 random set pairs for each of six calibrated size combinations. Sampling was uniform over the unique ChEMBL molecules, without replacement within each set and independently across sets; cross-set overlap remained allowed. The original model was held fixed. Exact zero-score P = 1 was retained.

At sizes 10 × 100, nominal P ≤ 0.01 occurred in 2.014% of independent null samples, and nominal P ≤ 0.001 occurred in 0.342%. The corresponding Bonferroni-adjusted excess-test P-values were approximately 7.68e-15 and 1.11e-6. There were 24 predeclared size/level checks. Other small-set combinations also showed excess rejection. A good pooled KS statistic therefore did not justify the deployed analytic approximation.

The new reference disables Gumbel extrapolation. Positive scores use the existing fixed-budget corrected Monte Carlo estimate `(k+1)/(B+1)` at the actual set sizes, with inclusive `>=` comparisons; zero scores retain exact P = 1. Existing screening/Balanced/Research budgets remain 999 / 199,999 / 999,999. Z is descriptive and is calculated from the actual size-specific sample mean and standard deviation; it is unavailable when variance is zero. Compatibility coefficient slots are placeholders (zero mean/standard-deviation amplitudes) and are not fitted MDDR/ChEMBL coefficients.

A separate empirical-null audit used 8,191 reference samples and 8,192 independent observed samples for each size pair. No excess rejection was detected across the 24 declared checks. Beta-binomial predictive tests account for uncertainty in the shared finite reference CDF; ties remain conservative. This is a check of the implemented uniform-set null, not a validation against biological inactive compounds.

## Leakage controls and retrieval comparison

Achiral Murcko scaffolds were hashed once for every source molecule. Hash buckets 2–4 supplied training references (264,076 molecules); bucket 1 supplied selection compounds (90,934), and bucket 0 supplied test compounds (86,850). All molecules sharing a scaffold, including the empty scaffold, remained in one bucket across every target. No held-out query scaffold remained in any training target set or the training background pool.

After retaining training target sets with at least 10 ligands, 1,299 targets remained for evaluation. Selection and test each used 512 distinct scaffolds. Target anchors were shuffled before choosing compounds, reducing domination by the largest ligand sets. Truth included every known active target for each query among the retained targets. Unmeasured pairs were never labelled inactive.

The 30 configurations were RDKit-path and Morgan radius 2 at thresholds 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.56, 0.57, 0.60, 0.65, 0.70, 0.75 and 0.80. Initial selection used 8,191 null samples. The winner was frozen before any test-query scoring. Confirmation used 32,767 samples, followed by the predeclared 199,999-sample numerical check. Study streams share prefixes across target sizes for efficiency; each marginal follows the same uniform-set null, but these streams are not bitwise identical to production's size-keyed streams.

MRR is the mean reciprocal rank of the first known active target. Recall@K is the mean fraction of each query's known active targets recovered among K results; Recall@1 is not top-1 accuracy. Paired bootstrap intervals resample 512 distinct query scaffolds 10,000 times. Shared targets can still introduce dependence, and these intervals do not cover dataset/assay-selection uncertainty.

| Null trials | Configuration | MRR | Recall@1 | Recall@5 | Recall@10 |
| ---: | --- | ---: | ---: | ---: | ---: |
| 32,767 | Existing path / 0.56 | 0.5566 | 0.3589 | 0.5596 | 0.6188 |
| 32,767 | Selected Morgan / 0.45 | 0.8344 | 0.5719 | 0.7935 | 0.8185 |
| 199,999 | Existing path / 0.56 | 0.5543 | 0.3569 | 0.5602 | 0.6167 |
| 199,999 | Selected Morgan / 0.45 | 0.8376 | 0.5761 | 0.7935 | 0.8185 |

At 199,999 trials, the paired MRR gain was **0.2833 (95% CI 0.2457–0.3228)** and Recall@10 gain was **0.2018 (0.1653–0.2398)**. Both predeclared adoption gates passed. The full 30-candidate selection, per-query results, audit counts and protocols are in [the immutable evidence report](chembl37-parameter-validation-results.json).

| Maximum-Tanimoto baseline | MRR | Recall@10 |
| --- | ---: | ---: |
| Morgan | 0.8798 | 0.8659 |
| RDKit-path | 0.8451 | 0.8437 |

A descriptive anchor-family breakdown appears in [family results](chembl37-family-validation-results.json). It was computed after the global decision and never used to tune parameters. Intervals are pointwise, and strata with fewer than 10 queries have no reported bootstrap interval.

| Anchor family | Queries | Existing MRR | Selected MRR |
| --- | ---: | ---: | ---: |
| GPCR | 74 | 0.5484 | 0.8462 |
| Kinase | 130 | 0.4956 | 0.8328 |
| Ion channel | 30 | 0.6943 | 0.8613 |
| Protease | 38 | 0.4733 | 0.8160 |
| Cytochrome P450 | 8 | 0.7505 | 0.6256 |
| Transporter | 15 | 0.8722 | 0.8671 |

**Improvement is not uniform across families:** Cytochrome P450 decreased on eight queries, and Transporter decreased slightly on 15. These small, post-hoc strata do not justify family-specific parameter choices; additional independent queries are needed.

The existing maximum-Tanimoto baseline is also reported. On this single-molecule holdout it outperformed the selected SEA configuration. Improving parameters within SEA does not establish that SEA is the best method for single-compound retrieval. Ensemble biological retrieval, measured inactive discrimination, rare families, alternative potency filters, chirality-aware representations and external-release validation remain separate work.

## Adoption and reproducibility

- Adopted dataset: `9d52be7ff211a45d2709`.
- Model: `59cb366620706b4e3870`.
- Evidence SHA256: `e8202ff920d88eb11c1bd3d498d861d41f433b32976b30b8f2fcb037f21d3f7c`.
- [Published manifest](chembl37-adopted-model.json) and [exact-cohort recipe](../config/chembl37-validated-parameters.json).

The new version retains the same 1,517 target sets, 441,860 standardized structures, activity filters, and source provenance. Only representation and statistical parameters change. The full reference is restored for application use; the scaffold-reduced reference exists only for evaluation. Old dataset IDs and job artifacts remain unchanged. The UI and preload selection prefer the approved version within ChEMBL 37, while old versions remain selectable. API/Python callers should explicitly pin a dataset ID.

Reproduction uses the locked Pixi environment and local source dataset:

```bash
pixi run python scripts/prepare_parameter_validation.py \
  --dataset-id e898df53d9bb5240f087 \
  --cache .cache/chembl37-parameter-validation --workers 4
cp config/chembl37-validation-protocol.json .cache/chembl37-parameter-validation/protocol.json
cp config/chembl37-production-protocol.json .cache/chembl37-parameter-validation/production-confirmation-protocol.json
pixi run python scripts/validate_chembl_parameters.py \
  --dataset-id e898df53d9bb5240f087 --cache .cache/chembl37-parameter-validation \
  --output .cache/chembl37-parameter-validation/retrieval-results.json
pixi run python scripts/confirm_parameter_budget.py \
  --dataset-id e898df53d9bb5240f087 --cache .cache/chembl37-parameter-validation \
  --output .cache/chembl37-parameter-validation/production-budget-results.json
pixi run python scripts/audit_chembl_null.py \
  --dataset-id e898df53d9bb5240f087 --cache .cache/chembl37-parameter-validation/null-audit \
  --output .cache/chembl37-parameter-validation/null-audit-results.json
pixi run python scripts/audit_empirical_null.py \
  --dataset-id e898df53d9bb5240f087 --cache .cache/chembl37-parameter-validation \
  --fingerprint Morgan --threshold 0.45 \
  --output .cache/chembl37-parameter-validation/empirical-null-results.json
```

The last command is specific to the frozen winner in this executed study; use the independently frozen parameters when reproducing on different data. Preparation resumes completed chunks. Larger null distributions require CPU time and several GB of temporary space. They create no application jobs or API keys. A rerun's timing metadata changes report hashes and can produce a different content-addressed model ID even when selected numerical parameters agree.

To adopt only after all criteria pass, stop existing services first so an older running worker cannot combine a new reference representation with its old query-fingerprint code:

```bash
./scripts/stop.sh
pixi run python scripts/adopt_chembl_parameters.py \
  --dataset-id e898df53d9bb5240f087 --cache .cache/chembl37-parameter-validation \
  --report docs/chembl37-parameter-validation-results.json
./scripts/start.sh
```

The deployment for this study was first created and checked in an isolated staging directory. All 30 completed result hashes, six API-key records and five original reference/metadata files were preserved when the idle lab services were restarted. Optional warming uses the adopted ID with `sea warm DATASET_ID --profile balanced --query-size 1 --limit 20`. Changing the sampling budget does not require rebuilding reference structures. Do not run the old KS-only `sea calibrate` command expecting it to preserve this validation decision: it creates a new experimental model, retains the reference fingerprint, and removes the approved label.

## Reuse the published decision for the exact cohort

If the original dataset `e898df53d9bb5240f087` is already installed, the published recipe rebuilds the verified reference without repeating the entire study:

```bash
./scripts/stop.sh
pixi run sea apply-parameters e898df53d9bb5240f087 \
  --recipe config/chembl37-validated-parameters.json
./scripts/start.sh
```

The command checks source identity, evidence checksum, adoption gates, the digest of **all** regenerated fingerprints and the expected new dataset identity before publishing. It clones existing annotations. A different cohort, filter set, source file or toolkit can require a new validation study; the recipe refuses mismatches. Standard `sea build` remains the original experimental builder and does not automatically certify another dataset or ChEMBL release.

## Interactive performance

An isolated API/two-worker run with this new reference and an empty null cache returned evidence in **0.707 s** and screening in **5.989 s**. A repeated search while refinement ran returned evidence in **0.401 s** and screening in **0.410 s**. Reference-loading readiness took another 6.053 s at startup. [Raw timings](chembl37-validated-performance-results.json) describe one local run; they are not percentile or matched-repeat speedup estimates. Full cold refinement was deliberately paused, so no final-refinement latency is claimed.

Outside interactive searches, the 20 most frequent target sizes were precomputed at 199,999 trials in staging for deployment. These caches are tied to the new immutable dataset ID and do not change its statistical model.

## Primary sources

- Keiser et al. (2007), [Relating protein pharmacology by ligand chemistry](https://doi.org/10.1038/nbt1284), and its Supplementary Methods: original MDDR/Daylight SEA framework.
- [RDKit fingerprint documentation](https://www.rdkit.org/docs/GettingStartedInPython.html#morgan-fingerprints-circular-fingerprints): explicit representation choices; implementation and toolkit version are recorded in the reference.
- Baldi and Nasr (2010), [When is Chemical Similarity Significant?](https://pmc.ncbi.nlm.nih.gov/articles/PMC2914517/): representation- and query-dependent similarity statistics. The model used here is the specified empirical uniform-set null, not that paper's analytic model.
