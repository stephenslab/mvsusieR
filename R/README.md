## Data

- `regression_data.R` implements several classes for different input data types and relevant algebraic operations on them.

## Model

Implementation of SuSiE model falls roughly in the structure of the SuSiE manuscript. 
That is, we introduce Bayesian regression model,
followed by Single Effect Regression (SER), and finally the SuSiE model.
Implementation-wise,

1. `*_regression.R` implements several classes for Bayesian regression models using different priors.
    - A conventional univariate multiple regression with normal prior.
    - A multivariate regression with fixed multivariate normal prior (prototype not optimized for production).
    - A multivariate regression with fixed multivariate normal mixture prior (the MASH model).
2. `single_effect_model.R` implements the SER model. It inherits Bayesian regression models.
3. `ibss_algorithm.R` implements the IBSS algorithm for SuSiE model.

## Misc

- `msusie.R` implements interface `msusie()` function.
- `utils.R` contains some utility functions.

## A note on `R6` class

`R6` is pretty [easy to learn](https://r6.r-lib.org/articles/Introduction.html) just by the length of its documentation.
However the constraint that [`private` cannot have same name as `public` and `active`](https://github.com/r-lib/R6/issues/200) is annoying.

Two quick pointers on `R6` class if you dont want to bother reading its documentation:

1. Convention `private$<name>` refers to private member (ie, variable) or method (ie, function); convention `self$<name>` refers to public member or method, or, active bindings. Please be careful which to use.
2. If you've got a class `A` to set `B` via `B=A` then modify `B`, you will also be modifying `A`! To prevent this from happening you've got to use this: `B=A$clone(deep=TRUE)`.

I use the following convention in my code:

1. There is no `public` member variable. All variables are private.
2. For variables meant to be used outside class I expose them to `active`.
    - By default active bindings for private variables are read-only, appear as `function()`; an attempt to assign values to them will raise an error.
    - Some active bindings will allow setting values of variables, and appear as `function(v)`.
    - Private members to be exposed to `active`, for the time being, will follow the convention of `.<name>` and the corresponding exposed variable will be `<name>`. This is the annoying constraint I was referring to above. Hopefully the R6 developers can change the behavior or convince me to believe it is good idea.

Otherwise, the only downside as I found using `R6` is that the trace back for errors is not as clear as without using it. But R language in general doesn't have good error trace back (compared to Python which can point to you the exact lines of codes that are problematic), so `R6` is not making it that much worse.