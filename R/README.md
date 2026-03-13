# mvsusieR R/ directory

File organization follows [susieR](https://github.com/stephenslab/susieR)
conventions so the two packages feel like parts of the same codebase.

## Entry points

- `mvsusie.R` -- `mvsusie()`, `mvsusie_rss()`, `mvsusie_ss()` 

## Data & constructors

- `mvsusie_constructors.R` -- `create_mvsusie_data`,
  `set_mvsusie_residual_variance`, `create_mvsusie_ss_data`, and
  model/fitted initialization functions.

## S3 methods (dispatched by susieR generics)

- `individual_data_methods.R` -- all `.mv_individual` methods.
- `sufficient_stats_methods.R` -- all `.mv_ss` methods.
- `model_methods.R` -- all `.mvsusie` model accessors.

## Utilities

- `mvsusie_utils.R` -- matrix helpers, numerical utilities, benchmarking tools,
  softmax, lfsr computation.

## Data

- `example_dataset.R` -- `mvsusie_sim1` and `simdata` documentation.

## Infrastructure

- `zzz.R` -- `.onLoad` for S3 method registration and susieR / mashr function caching.
