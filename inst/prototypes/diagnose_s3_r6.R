#!/usr/bin/env Rscript
# =============================================================================
# Diagnose S3 vs R6 discrepancies for V4, V8a, V8b
# =============================================================================
# Strategy: run S3 with R6 defaults to isolate default vs code difference

cat("=== Diagnosing S3 vs R6 Discrepancies ===\n\n")

suppressMessages(library(pkgload))

# Load R6 reference from local patched copy (check_null_threshold warmup disabled)
ref_source <- normalizePath("~/GIT/R6_mvsusieR")
ref_env <- suppressMessages(pkgload::load_all(ref_source, export_all = TRUE,
                                               quiet = TRUE))

# Load S3 development package
dev_source <- normalizePath(file.path(getwd()))
dev_env <- suppressMessages(pkgload::load_all(dev_source, export_all = TRUE,
                                               quiet = TRUE))

# R6 mvsusie helper
r6_mvsusie <- function(...) {
  args <- list(...)
  if ("precompute_cache" %in% names(args)) {
    args$precompute_covariances <- args$precompute_cache
    args$precompute_cache <- NULL
  }
  args$estimate_prior_mixture_weights <- NULL
  args$mixture_weight_method <- NULL
  do.call(ref_env$env$mvsusie, args)
}

r6_create_mixture_prior <- function(...) ref_env$env$create_mixture_prior(...)

# =============================================================================
# V4: linkage_vs_pleiotropy — DIAGNOSIS
# =============================================================================
cat("\n===================================================================\n")
cat("V4: LINKAGE vs PLEIOTROPY DIAGNOSIS\n")
cat("===================================================================\n\n")

set.seed(1)
data(N3finemapping, package = "susieR")
X <- N3finemapping$X
n <- nrow(X); p <- ncol(X); r <- 10
R_mat <- cor(X)
pairs <- which(abs(R_mat) > 0.65 & abs(R_mat) < 0.75 & upper.tri(R_mat),
               arr.ind = TRUE)
snp1 <- pairs[1, 1]; snp2 <- pairs[1, 2]
cat(sprintf("Causal SNPs: %d, %d (LD = %.3f)\n", snp1, snp2,
            R_mat[snp1, snp2]))
B4 <- matrix(0, p, r)
B4[snp1, 1:5] <- rnorm(5, 0, 1)
B4[snp2, 6:10] <- rnorm(5, 0, 1)
Y4 <- X %*% B4 + matrix(rnorm(n * r), n, r)

prior_s3 <- dev_env$env$create_mixture_prior(R = r)
prior_r6 <- r6_create_mixture_prior(R = r)

# 1. S3 with S3 defaults
cat("\n--- 1. S3 with S3 defaults ---\n")
s3_default <- dev_env$env$mvsusie(X, Y4, L = 10, prior_variance = prior_s3,
                                   estimate_prior_variance = TRUE,
                                   standardize = TRUE, tol = 0.01)
cat(sprintf("  CS: %d, PIP@snp1: %.4f, PIP@snp2: %.4f\n",
            length(s3_default$sets$cs),
            s3_default$pip[snp1], s3_default$pip[snp2]))
cat(sprintf("  Coef vs truth: %.4f\n",
            cor(as.vector(B4), as.vector(coef(s3_default)[-1, ]))))
cat(sprintf("  ELBO: %.2f, niter: %d\n",
            tail(s3_default$elbo, 1), s3_default$niter))

# 2. R6 with R6 defaults
cat("\n--- 2. R6 with R6 defaults ---\n")
r6_default <- r6_mvsusie(X, Y4, L = 10, prior_variance = prior_r6,
                          estimate_prior_variance = TRUE,
                          standardize = TRUE, tol = 0.01,
                          precompute_cache = FALSE,
                          estimate_prior_method = "EM",
                          estimate_residual_variance = FALSE)
cat(sprintf("  CS: %d, PIP@snp1: %.4f, PIP@snp2: %.4f\n",
            length(r6_default$sets$cs),
            r6_default$pip[snp1], r6_default$pip[snp2]))
cat(sprintf("  Coef vs truth: %.4f\n",
            cor(as.vector(B4), as.vector(r6_default$coef[-1, ]))))
cat(sprintf("  ELBO: %.2f, niter: %d\n",
            tail(r6_default$elbo, 1), r6_default$niter))

# 3. S3 with R6 defaults (isolate: is it defaults or code?)
cat("\n--- 3. S3 with R6 defaults ---\n")
s3_r6def <- dev_env$env$mvsusie(X, Y4, L = 10, prior_variance = prior_s3,
                                 estimate_prior_variance = TRUE,
                                 estimate_prior_method = "EM",
                                 estimate_residual_variance = FALSE,
                                 precompute_cache = FALSE,
                                 standardize = TRUE, tol = 0.01)
cat(sprintf("  CS: %d, PIP@snp1: %.4f, PIP@snp2: %.4f\n",
            length(s3_r6def$sets$cs),
            s3_r6def$pip[snp1], s3_r6def$pip[snp2]))
cat(sprintf("  PIP cor to R6: %.6f\n", cor(s3_r6def$pip, r6_default$pip)))
cat(sprintf("  Coef vs truth: %.4f\n",
            cor(as.vector(B4), as.vector(coef(s3_r6def)[-1, ]))))
cat(sprintf("  ELBO: %.2f, niter: %d\n",
            tail(s3_r6def$elbo, 1), s3_r6def$niter))

# 4. Isolate each default: S3 with only estimate_residual_variance=FALSE
cat("\n--- 4. S3 with estimate_residual_variance=FALSE (other S3 defaults) ---\n")
s3_noerv <- dev_env$env$mvsusie(X, Y4, L = 10, prior_variance = prior_s3,
                                 estimate_prior_variance = TRUE,
                                 estimate_residual_variance = FALSE,
                                 standardize = TRUE, tol = 0.01)
cat(sprintf("  CS: %d, PIP@snp1: %.4f, PIP@snp2: %.4f\n",
            length(s3_noerv$sets$cs),
            s3_noerv$pip[snp1], s3_noerv$pip[snp2]))
cat(sprintf("  PIP cor to R6: %.6f\n", cor(s3_noerv$pip, r6_default$pip)))
cat(sprintf("  Coef vs truth: %.4f\n",
            cor(as.vector(B4), as.vector(coef(s3_noerv)[-1, ]))))

# 5. S3 with only estimate_prior_method="EM"
cat("\n--- 5. S3 with estimate_prior_method='EM' (other S3 defaults) ---\n")
s3_em <- dev_env$env$mvsusie(X, Y4, L = 10, prior_variance = prior_s3,
                              estimate_prior_variance = TRUE,
                              estimate_prior_method = "EM",
                              standardize = TRUE, tol = 0.01)
cat(sprintf("  CS: %d, PIP@snp1: %.4f, PIP@snp2: %.4f\n",
            length(s3_em$sets$cs),
            s3_em$pip[snp1], s3_em$pip[snp2]))
cat(sprintf("  PIP cor to S3 default: %.6f\n", cor(s3_em$pip, s3_default$pip)))
cat(sprintf("  Coef vs truth: %.4f\n",
            cor(as.vector(B4), as.vector(coef(s3_em)[-1, ]))))

cat("\n--- V4 DIAGNOSIS ---\n")
cat("Which default change matters most?\n")
cat(sprintf("  S3 defaults:       coef_vs_truth = %.4f\n",
            cor(as.vector(B4), as.vector(coef(s3_default)[-1, ]))))
cat(sprintf("  S3 + ERV=FALSE:    coef_vs_truth = %.4f\n",
            cor(as.vector(B4), as.vector(coef(s3_noerv)[-1, ]))))
cat(sprintf("  S3 + EM:           coef_vs_truth = %.4f\n",
            cor(as.vector(B4), as.vector(coef(s3_em)[-1, ]))))
cat(sprintf("  S3 + R6 defaults:  coef_vs_truth = %.4f\n",
            cor(as.vector(B4), as.vector(coef(s3_r6def)[-1, ]))))
cat(sprintf("  R6 defaults:       coef_vs_truth = %.4f\n",
            cor(as.vector(B4), as.vector(r6_default$coef[-1, ]))))


# =============================================================================
# V8a: UV-MV agreement — DIAGNOSIS
# =============================================================================
cat("\n\n===================================================================\n")
cat("V8a: UV-MV AGREEMENT DIAGNOSIS\n")
cat("===================================================================\n\n")

set.seed(1)
sim8 <- dev_env$env$mvsusie_sim1(r = 10)
r <- ncol(sim8$y)
causal8 <- which(rowSums(abs(sim8$b)) > 0)
U_diag <- 0.2 * diag(r)
prior_s3 <- dev_env$env$create_mixture_prior(
  list(matrices = list(U_diag), weights = 1))
prior_r6 <- r6_create_mixture_prior(
  list(matrices = list(U_diag), weights = 1))

cat(sprintf("True causal SNPs: %s\n", paste(causal8, collapse = ", ")))

# 1. S3 with S3 defaults
cat("\n--- 1. S3 with S3 defaults ---\n")
s3_8d <- dev_env$env$mvsusie(sim8$X, sim8$y, L = 10, prior_variance = prior_s3,
                              estimate_prior_variance = TRUE,
                              standardize = TRUE, tol = 0.01)
cat(sprintf("  CS: %d, causal covered: %d/%d\n",
            length(s3_8d$sets$cs),
            sum(sapply(causal8, function(j) any(sapply(s3_8d$sets$cs,
                                                       function(cs) j %in% cs)))),
            length(causal8)))
cat(sprintf("  Coef vs truth: %.4f\n",
            cor(as.vector(sim8$b), as.vector(coef(s3_8d)[-1, ]))))
cat(sprintf("  ELBO: %.2f, niter: %d\n",
            tail(s3_8d$elbo, 1), s3_8d$niter))
cat(sprintf("  CS sizes: %s\n",
            paste(sapply(s3_8d$sets$cs, length), collapse = ", ")))

# 2. R6 with R6 defaults
cat("\n--- 2. R6 with R6 defaults ---\n")
r6_8d <- r6_mvsusie(sim8$X, sim8$y, L = 10, prior_variance = prior_r6,
                     estimate_prior_variance = TRUE,
                     standardize = TRUE, tol = 0.01,
                     precompute_cache = FALSE,
                     estimate_prior_method = "EM",
                     estimate_residual_variance = FALSE)
cat(sprintf("  CS: %d, causal covered: %d/%d\n",
            length(r6_8d$sets$cs),
            sum(sapply(causal8, function(j) any(sapply(r6_8d$sets$cs,
                                                       function(cs) j %in% cs)))),
            length(causal8)))
cat(sprintf("  Coef vs truth: %.4f\n",
            cor(as.vector(sim8$b), as.vector(r6_8d$coef[-1, ]))))
cat(sprintf("  ELBO: %.2f, niter: %d\n",
            tail(r6_8d$elbo, 1), r6_8d$niter))
cat(sprintf("  CS sizes: %s\n",
            paste(sapply(r6_8d$sets$cs, length), collapse = ", ")))

# 3. S3 with R6 defaults
cat("\n--- 3. S3 with R6 defaults ---\n")
s3_8r6 <- dev_env$env$mvsusie(sim8$X, sim8$y, L = 10, prior_variance = prior_s3,
                                estimate_prior_variance = TRUE,
                                estimate_prior_method = "EM",
                                estimate_residual_variance = FALSE,
                                precompute_cache = FALSE,
                                standardize = TRUE, tol = 0.01)
cat(sprintf("  CS: %d, PIP cor to R6: %.6f\n",
            length(s3_8r6$sets$cs), cor(s3_8r6$pip, r6_8d$pip)))
cat(sprintf("  Coef vs truth: %.4f\n",
            cor(as.vector(sim8$b), as.vector(coef(s3_8r6)[-1, ]))))
cat(sprintf("  ELBO: %.2f, niter: %d\n",
            tail(s3_8r6$elbo, 1), s3_8r6$niter))
cat(sprintf("  CS sizes: %s\n",
            paste(sapply(s3_8r6$sets$cs, length), collapse = ", ")))

# 4. S3 with estimate_residual_variance=FALSE only
cat("\n--- 4. S3 with ERV=FALSE only ---\n")
s3_8noerv <- dev_env$env$mvsusie(sim8$X, sim8$y, L = 10, prior_variance = prior_s3,
                                  estimate_prior_variance = TRUE,
                                  estimate_residual_variance = FALSE,
                                  standardize = TRUE, tol = 0.01)
cat(sprintf("  CS: %d, PIP cor to R6: %.6f\n",
            length(s3_8noerv$sets$cs), cor(s3_8noerv$pip, r6_8d$pip)))
cat(sprintf("  Coef vs truth: %.4f\n",
            cor(as.vector(sim8$b), as.vector(coef(s3_8noerv)[-1, ]))))

# 5. S3 with EM only
cat("\n--- 5. S3 with EM only ---\n")
s3_8em <- dev_env$env$mvsusie(sim8$X, sim8$y, L = 10, prior_variance = prior_s3,
                               estimate_prior_variance = TRUE,
                               estimate_prior_method = "EM",
                               standardize = TRUE, tol = 0.01)
cat(sprintf("  CS: %d, PIP cor to S3 default: %.6f\n",
            length(s3_8em$sets$cs), cor(s3_8em$pip, s3_8d$pip)))
cat(sprintf("  Coef vs truth: %.4f\n",
            cor(as.vector(sim8$b), as.vector(coef(s3_8em)[-1, ]))))

cat("\n--- V8a DIAGNOSIS ---\n")
cat(sprintf("  S3 defaults:       CS=%d, coef_vs_truth = %.4f\n",
            length(s3_8d$sets$cs),
            cor(as.vector(sim8$b), as.vector(coef(s3_8d)[-1, ]))))
cat(sprintf("  S3 + ERV=FALSE:    CS=%d, coef_vs_truth = %.4f\n",
            length(s3_8noerv$sets$cs),
            cor(as.vector(sim8$b), as.vector(coef(s3_8noerv)[-1, ]))))
cat(sprintf("  S3 + EM:           CS=%d, coef_vs_truth = %.4f\n",
            length(s3_8em$sets$cs),
            cor(as.vector(sim8$b), as.vector(coef(s3_8em)[-1, ]))))
cat(sprintf("  S3 + R6 defaults:  CS=%d, coef_vs_truth = %.4f\n",
            length(s3_8r6$sets$cs),
            cor(as.vector(sim8$b), as.vector(coef(s3_8r6)[-1, ]))))
cat(sprintf("  R6 defaults:       CS=%d, coef_vs_truth = %.4f\n",
            length(r6_8d$sets$cs),
            cor(as.vector(sim8$b), as.vector(r6_8d$coef[-1, ]))))


# =============================================================================
# V8b: ELBO/convergence — DIAGNOSIS
# =============================================================================
cat("\n\n===================================================================\n")
cat("V8b: ELBO/CONVERGENCE DIAGNOSIS\n")
cat("===================================================================\n\n")

prior_s3_canon <- dev_env$env$create_mixture_prior(R = 10)
prior_r6_canon <- r6_create_mixture_prior(R = 10)

# 1. S3 with S3 defaults
cat("--- 1. S3 with S3 defaults ---\n")
s3_8bd <- dev_env$env$mvsusie(sim8$X, sim8$y, L = 10,
                               prior_variance = prior_s3_canon,
                               estimate_prior_variance = TRUE,
                               standardize = TRUE, tol = 1e-3)
cat(sprintf("  CS: %d, coef_vs_truth: %.4f, niter: %d\n",
            length(s3_8bd$sets$cs),
            cor(as.vector(sim8$b), as.vector(coef(s3_8bd)[-1, ])),
            s3_8bd$niter))

# 2. R6 defaults
cat("--- 2. R6 with R6 defaults ---\n")
r6_8bd <- r6_mvsusie(sim8$X, sim8$y, L = 10, prior_variance = prior_r6_canon,
                      estimate_prior_variance = TRUE,
                      standardize = TRUE, tol = 1e-3,
                      precompute_cache = FALSE,
                      estimate_prior_method = "EM",
                      estimate_residual_variance = FALSE)
cat(sprintf("  CS: %d, coef_vs_truth: %.4f, niter: %d\n",
            length(r6_8bd$sets$cs),
            cor(as.vector(sim8$b), as.vector(r6_8bd$coef[-1, ])),
            r6_8bd$niter))

# 3. S3 with R6 defaults
cat("--- 3. S3 with R6 defaults ---\n")
s3_8br6 <- dev_env$env$mvsusie(sim8$X, sim8$y, L = 10,
                                 prior_variance = prior_s3_canon,
                                 estimate_prior_variance = TRUE,
                                 estimate_prior_method = "EM",
                                 estimate_residual_variance = FALSE,
                                 precompute_cache = FALSE,
                                 standardize = TRUE, tol = 1e-3)
cat(sprintf("  CS: %d, PIP cor to R6: %.6f, coef_vs_truth: %.4f, niter: %d\n",
            length(s3_8br6$sets$cs),
            cor(s3_8br6$pip, r6_8bd$pip),
            cor(as.vector(sim8$b), as.vector(coef(s3_8br6)[-1, ])),
            s3_8br6$niter))

# 4. S3 with ERV=FALSE only
cat("--- 4. S3 with ERV=FALSE only ---\n")
s3_8bnoerv <- dev_env$env$mvsusie(sim8$X, sim8$y, L = 10,
                                    prior_variance = prior_s3_canon,
                                    estimate_prior_variance = TRUE,
                                    estimate_residual_variance = FALSE,
                                    standardize = TRUE, tol = 1e-3)
cat(sprintf("  CS: %d, coef_vs_truth: %.4f\n",
            length(s3_8bnoerv$sets$cs),
            cor(as.vector(sim8$b), as.vector(coef(s3_8bnoerv)[-1, ]))))


# =============================================================================
# FINAL DIAGNOSIS SUMMARY
# =============================================================================
cat("\n\n===================================================================\n")
cat("FINAL DIAGNOSIS SUMMARY\n")
cat("===================================================================\n\n")

cat("V4 linkage_vs_pleiotropy:\n")
cat(sprintf("  Is it a defaults issue? S3+R6_defaults PIP cor to R6 = %.6f\n",
            cor(s3_r6def$pip, r6_default$pip)))
v4_defaults <- cor(s3_r6def$pip, r6_default$pip) > 0.99
cat(sprintf("  Conclusion: %s\n",
            ifelse(v4_defaults,
                   "YES - difference is entirely from defaults, not a code bug",
                   "NO - there may be a code difference")))
cat(sprintf("  Which default matters? ERV=FALSE coef_truth=%.4f, EM coef_truth=%.4f, both=%.4f\n",
            cor(as.vector(B4), as.vector(coef(s3_noerv)[-1, ])),
            cor(as.vector(B4), as.vector(coef(s3_em)[-1, ])),
            cor(as.vector(B4), as.vector(coef(s3_r6def)[-1, ]))))

cat("\nV8a UV-MV agreement:\n")
cat(sprintf("  Is it a defaults issue? S3+R6_defaults PIP cor to R6 = %.6f\n",
            cor(s3_8r6$pip, r6_8d$pip)))
v8a_defaults <- cor(s3_8r6$pip, r6_8d$pip) > 0.99
cat(sprintf("  Conclusion: %s\n",
            ifelse(v8a_defaults,
                   "YES - difference is entirely from defaults, not a code bug",
                   "NO - there may be a code difference")))
cat(sprintf("  S3 defaults CS=%d vs R6 defaults CS=%d\n",
            length(s3_8d$sets$cs), length(r6_8d$sets$cs)))
cat(sprintf("  S3+R6_defaults CS=%d\n", length(s3_8r6$sets$cs)))
cat(sprintf("  Key: ERV=FALSE → CS=%d (matches R6? %s)\n",
            length(s3_8noerv$sets$cs),
            ifelse(length(s3_8noerv$sets$cs) == length(r6_8d$sets$cs), "YES", "NO")))

cat("\nV8b ELBO/convergence:\n")
cat(sprintf("  Is it a defaults issue? S3+R6_defaults PIP cor to R6 = %.6f\n",
            cor(s3_8br6$pip, r6_8bd$pip)))
v8b_defaults <- cor(s3_8br6$pip, r6_8bd$pip) > 0.99
cat(sprintf("  Conclusion: %s\n",
            ifelse(v8b_defaults,
                   "YES - difference is entirely from defaults, not a code bug",
                   "NO - there may be a code difference")))
cat(sprintf("  S3 defaults CS=%d vs R6 defaults CS=%d\n",
            length(s3_8bd$sets$cs), length(r6_8bd$sets$cs)))
cat(sprintf("  S3+R6_defaults CS=%d\n", length(s3_8br6$sets$cs)))

cat("\n\nOverall:\n")
all_defaults <- v4_defaults && v8a_defaults && v8b_defaults
if (all_defaults) {
  cat("ALL discrepancies are from default parameter differences, NOT code bugs.\n")
  cat("The S3 refactoring is correct — only defaults differ.\n\n")
  cat("Default differences and their impact:\n")
  cat("  estimate_residual_variance: S3=TRUE (re-estimates V) vs R6=FALSE (fixed V)\n")
  cat("    → Biggest impact on pruning behavior. With ERV=TRUE, more effects survive.\n")
  cat("  estimate_prior_method: S3='optim' vs R6='EM'\n")
  cat("    → Usually minor impact. optim is faster but may converge differently.\n")
} else {
  cat("Some discrepancies may involve code differences — further investigation needed.\n")
}
