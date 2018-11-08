## Data

- `regression_data.R` implements several classes for different input data types and relevant algebraic operations on them.

## Model

Implementation of SuSiE model falls roughly in the structure of the SuSiE manuscript. 
That is, we introduce Bayesian multiple regression model (BMR), 
followed by Single Effect BMR (SER), and finally the SuSiE model.
Implementation-wise,

1. `bayesian_multiple_regression.R` implements several classes for BMR models using different priors.
    - A conventional univariate multiple regression with normal prior
    - A multivariate regression with fixed multivariate normal prior (the MASH model)
2. `single_effect_regression.R` implements the SER model. It inherits BMR models.
3. `susie_regression.R` implements the SuSiE model and algorithm.

## Misc

- `main.R` implements interface.
- `utils.R` contains some utility functions.