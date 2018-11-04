## Data

- `regression_data.R` implements several classes for different input data types and relevant algebraic operations on them.

## Model

Implementation of SuSiE model falls roughly in the structure of the SuSiE manuscript. 
That is, we introduce Bayesian multiple regression model (BMR), 
followed by Single Effect BMR (SER), and finally the SuSiE model.
Implementation-wise,

1. `bayesian_multiple_regression.R` implements several classes for Bayesian Multiple regression models using different priors.
2. `single_effect_regression.R` implements the SER model.
3. `susie_regression.R` implements the SuSiE model and algorithm.

Additionally we have a dedicated file for SuSiE CS related implementation,

- `susie_cs.R` implements computation and visualization of SuSiE credible sets.

## Misc

- `main.R` implements interface.
- `utils.R` contains some utility functions.