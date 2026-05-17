# Log: Fixing estimate_residual_variance for 3D missing data path

## Problem
When `estimate_residual_variance=TRUE` with missing data using the 3D path,
the ELBO is non-monotone because:
1. `model$bxxb` uses scalar `d[j] = N-1`, but in 3D path the per-outcome
   Gram matrix diagonal `G_j[r,r] = n_r - 1` (where n_r is per-outcome
   observation count). For 20% MCAR, this overestimates by ~25%.
2. After V is updated, stale V-dependent quantities (Vinv, elbo_const, svs)
   cause incorrect ELBO and SER computations.

## Injection fixes (Changes A-E) — DONE
Solved problem #2: inject updated V-dependent quantities from model into data
at each call site (compute_residuals, compute_ser_statistics, compute_elbo_3d,
estimate_residual_variance_mv, update_model_variance). All 1465 tests pass.

## Math background

### EM M-step for V with missing data
The M-step maximizes E_q[log p(Y_obs | B, V)] w.r.t. V:

log p(Y_obs | B, V) = sum_i [-|obs_i|/2 log(2pi) - 1/2 log|V[obs_i,obs_i]|
                      - 1/2 e_i' V[obs_i,obs_i]^{-1} e_i]

where e_i = Y_i[obs_i] - X_3d[i,:,obs_i]' B[:,obs_i].

Setting d/dV = 0 does NOT have closed form because V appears in sub-blocks
V[obs_k,obs_k] for different missingness patterns k.

### What the bxxb correction represents
For single effect l with gamma_l (selected variable):
  (X_3d b_l)_i[r] = X_3d[i,gamma_l,r] * beta_l[r]

So: [(X_3d b_l)' (X_3d b_l)]_{r,s} = beta_l[r] * beta_l[s] * G_{gamma_l}[r,s]

Taking expectation over q(gamma_l, beta_l):
  E_q[(X_3d b_l)' (X_3d b_l)]_{r,s} = sum_j alpha_lj * G_j[r,s] * mu2_lj[r,s]
This is a HADAMARD (element-wise) product of G_j and mu2_lj.

The posterior correction per effect is:
  correction_l = E_q[(X_3d b_l)' (X_3d b_l)] - (X_3d E[b_l])' (X_3d E[b_l])

## Attempt 1: Inner EM with Hadamard + matrix-mult bxxb (WRONG MATH)
- Used `tcrossprod(mu_l[j,], G_j %*% mu_l[j,])` = mu * mu' * G_j' (matrix product)
- Should be `G_j * tcrossprod(mu_l[j,])` (Hadamard product)
- Result: ELBO wildly divergent, V explodes to 27

## Attempt 2: Inner EM with correct Hadamard product
- Fixed to use `G_j * mu_outer` (element-wise *)
- Result: still wildly divergent; V oscillates between small and large values
- ELBO swings: +3570, -3290, +91700, -94700, ...

Analysis: The inner EM adds the posterior correction AFTER imputing missing
residuals. But the correction only covers observed entries (X_3d = 0 for missing).
The imputed residuals get amplified correction through Lambda_k, creating
instability.

## Attempt 3: Pairwise complete-case estimator (no inner EM)
V[r,s] = (crossprod(R_hat)[r,s] + correction[r,s]) / N_{rs}

where:
- R_hat = Y - X_3d E[B], missing entries = 0
- correction = sum_l (bxxb_3d_l - b1_XtX_b1_l) using Hadamard formula
- N_{rs} = number of samples where both r and s are observed

This is a method-of-moments estimator. It directly computes:
  V[r,s] = (1/N_{rs}) sum_{i: obs r,s} E_q[e_i[r] * e_i[s]]

Result: sigma2 values are REASONABLE (~1-2.5), but ELBO oscillates wildly
(+-30000+) when combined with estimate_prior_variance=TRUE.

## Attempt 4: Damped sigma2 update (eta=0.5)
- Added `model$sigma2 <- (1-eta)*model$sigma2 + eta*sigma2_candidate`
- Result: still diverges; damping delays but doesn't prevent instability.

## Attempt 5: Diagonal-exact + shrunk off-diagonal
- Diagonal from exact maximizer, off-diagonal from pairwise + 50% shrinkage
- Result: WORSE — V_prior explodes to 3.27e+06, ELBO becomes NaN,
  3006 Cholesky failures.

## Root cause analysis
- est_rv ALONE (est_pv=FALSE): ELBO monotone, converges in 14 iters
- est_pv ALONE (est_rv=FALSE): ELBO monotone, converges in 3 iters
- BOTH together: ELBO oscillates wildly

The problem is the INTERACTION between sigma2 and prior variance updates:
- sigma2 changes -> svs_inv changes NON-UNIFORMLY (per-pattern, per-variable)
- Non-uniform svs_inv perturbation causes prior variance optimizer to overshoot
- Overshooting causes ELBO swings

In the no-missing case, svs_inv = V^{-1} * d[j] changes UNIFORMLY (scalar
factor on all variables), so the optimizer adjusts smoothly. In the 3D path,
different missingness patterns create variable-specific svs_inv values that
change heterogeneously when sigma2 is updated.

## Attempt 6: Outer loop (SOLUTION) — DONE
Instead of updating sigma2 within the IBSS loop (where it destabilizes the
prior variance optimizer), use an outer loop:

1. Phase 1: Run IBSS with FIXED sigma2 until convergence (guaranteed monotone)
2. Phase 2: Compute sigma2 from converged posterior via pairwise estimator
3. Phase 3: Update data with new sigma2, fresh-start next IBSS run
4. Repeat until sigma2 converges (rel_change < tol) or max 10 outer iterations

Implementation in `mvsusie_core` (R/mvsusie.R):
- `use_outer_rv_loop <- estimate_residual_variance && !is.null(data$miss3d)`
- Inner loop runs `mvsusie_workhorse` with `estimate_residual_variance=FALSE`
- After convergence, `estimate_residual_variance_3d` computes new sigma2
- `set_residual_variance_3d` updates all V-dependent quantities in data
- No warm-starting (cleanup_model prunes mixture components, incompatible
  with original prior_variance on re-init)

Results:
- 20% MCAR: outer loop converges in 2 iterations, ELBO monotone within
  each inner loop, sigma2 diag = [1.035, 1.066, 1.073, 1.15, 1.005]
- 50% MCAR: outer loop converges in 2 iterations, sigma2 reasonable
- Block missing: outer loop converges in 3 iterations
- No missing: standard path unchanged, ELBO monotone
- PIP correlation between est_rv and no-est_rv: 0.9996
- All 1465 tests pass with 0 failures

### Files modified:
1. `R/mvsusie.R`: Added outer loop in `mvsusie_core` for 3D + est_rv case
2. `R/missing_y_utils.R`: Added `estimate_residual_variance_3d` function
   (pairwise complete-case estimator with Hadamard-product posterior correction)
3. `R/individual_data_methods.R`: `estimate_residual_variance_mv` dispatches
   to `estimate_residual_variance_3d` for 3D path; `update_model_variance`
   recomputes all V-dependent quantities after sigma2 update
