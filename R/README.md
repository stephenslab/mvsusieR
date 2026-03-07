# mvsusieR R/ directory

File organization follows [susieR](https://github.com/stephenslab/susieR)
conventions so the two packages feel like parts of the same codebase.

## Entry points

- `mvsusie.R` -- `mvsusie()`, `mvsusie_rss()`, `mvsusie_suff_stat()` and
  their internal S3 counterparts `mvsusie_s3()`, `mvsusie_suff_stat_s3()`.

## Data & constructors

- `mvsusie_constructors.R` -- `create_mvsusie_data`,
  `set_mvsusie_residual_variance`, `create_mvsusie_ss_data`, and
  model/fitted initialization functions.

## S3 methods (dispatched by susieR generics)

- `individual_data_methods.R` -- all `.mv_individual` methods.
- `sufficient_stats_methods.R` -- all `.mv_ss` methods.
- `model_methods.R` -- all `.mvsusie` model accessors.

## Algorithm

- `mvsusie_workhorse.R` -- the IBSS loop (calls `susieR::susie_workhorse`).

## Prior

- `mash_prior.R` -- `create_mash_prior` and prior helpers.
- `mash.R` -- `create_mixture_prior` (user-facing wrapper).

## Output & prediction

- `mvsusie_get_functions.R` -- `format_mvsusie_output`,
  `apply_mvsusie_dimnames`, `multivariate_lbf`.
- `predict.mvsusie.R` -- `coef.mvsusie`, `predict.mvsusie`.
- `mvsusie_plot.R` -- plotting functions.

## Utilities

- `mvsusie_utils.R` -- matrix helpers, numerical utilities, benchmarking tools,
  softmax, lfsr computation.

## Data

- `example_dataset.R` -- `mvsusie_sim1` and `simdata` documentation.

## Infrastructure

- `zzz.R` -- `.onLoad` for S3 method registration and mashr function caching.
