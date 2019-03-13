#' @title compute value_j * weight_j / sum(value_j * weight_j)
#' @keywords internal
safe_compute_weight = function(value, weight, log = TRUE) {
    if (!log) value = log(value)
    mvalue = max(value)
    w = exp(value-mvalue)
    w_weighted = w * weight
    weighted_sum_w = sum(w_weighted)
    return(list(alpha = as.vector(w_weighted / weighted_sum_w), log_total = log(weighted_sum_w) + mvalue))
}

#' @title SuSiE model extractor
#' @importFrom abind abind
#' @keywords internal
report_susie_model = function(d, m) {
    if (is.null(dim(m$posterior_b1[[1]]))) {
      mu = t(do.call(cbind, m$posterior_b1))
      mu2 = t(do.call(cbind, m$posterior_b2))
      b = colSums(t(m$pip)*mu)
    } else {
      mu =  aperm(abind::abind(m$posterior_b1,along=3), c(3,1,2))
      mu2 = aperm(abind::abind(m$posterior_b2,along=3), c(3,1,2))
      b = do.call(cbind, lapply(1:dim(mu)[3], function(i) colSums(t(m$pip) * mu[,,i])))
    }
    if (is.null(m$mixture_posterior_weights)) mixture_weights = NA
    else mixture_weights = aperm(abind::abind(m$mixture_posterior_weights,along=3), c(3,1,2))
    lfsr = aperm(abind::abind(m$lfsr,along=3), c(3,1,2))
    s = list(
        alpha = t(m$pip),
        mu = mu,
        mu2 = mu2,
        KL = m$kl,
        lbf = m$lbf,
        sigma2 = m$residual_variance,
        V = m$prior_variance,
        elbo = m$get_objective(dump=TRUE),
        niter = m$get_niter(),
        fitted = d$fitted,
        coef = d$rescale_coef(b),
        null_index = -9,
        mixture_weights = mixture_weights,
        lfsr = lfsr 
        )
    s$intercept = s$coef[1]
    class(s) = 'susie'
    return(s)
}

#' @title Compute condition specific posterior inclusion probability
#' @param m M&M model
#' @param prior_obj prior mixture object
#' @return P by R matrix of PIP per condition
#' @keywords internal
mmbr_get_pip_per_condition = function(m, prior_obj) {
  condition_pip = mmbr_get_alpha_per_condition(m, prior_obj)
  return(do.call(cbind, lapply(1:dim(condition_pip)[3], function(r) apply(condition_pip[,,r], 2, function(x) 1-prod(1-x)))))
}

#' @title Compute condition specific posterior inclusion probability per effect
#' @keywords internal
mmbr_get_alpha_per_condition = function(m, prior_obj) {
  condition_indicator = do.call(rbind, lapply(1:length(prior_obj$prior_covariance$xUlist), function(i) as.integer(diag(prior_obj$prior_covariance$xUlist[[i]]) != 0)))
  condition_pip = array(0, dim=dim(m$mu))
  for (r in 1:dim(condition_pip)[3]) {
    for (p in 1:length(condition_indicator[,r])) {
        condition_pip[,,r] = condition_pip[,,r] + m$mixture_weights[,,p] * condition_indicator[p,r]
    }
    condition_pip[,,r] = condition_pip[,,r] * m$alpha
  }
  return(condition_pip)
}

#' @title A null progressbar, because currently `progressbar_enabled` feature does not work for `progress_bar`
#' @importFrom R6 R6Class
#' @keywords internal
null_progress_bar = R6Class('null_progress_bar', public = list(tick = function(...) {}))

#' @title check if all elements are the same in matrix of J by R, J >> R
#' @keywords internal
is_mat_common = function(mat) {
  all((t(mat) - mat[1,]) == 0)
}

#' @title A simple simulation function to simulate some test data
#' @param n number of samples
#' @param p number of features
#' @param r number of conditions
#' @param s number of effect variables per condition if greater than 1; otherwise percentage of effect variables per condition
#' @param center_scale FALSE by default
#' @export
mmbr_sim1 = function(n=200,p=500,r=2,s=4,center_scale=FALSE) {
  X = matrix(rnorm(n*p,0,1),n,p)
  if (s>1) {
    beta = matrix(0, p, r)
    for (i in 1:r) beta[sample(1:p,s), i] = 1
  } else {
    beta = matrix(runif(p*r)>s, p, r)
  }
  y = X %*% beta + do.call(cbind, lapply(1:r, function(i) rnorm(n)))
  if (center_scale) {
    X = scale(X)
    y = y - apply(y,2,mean)
  }
  scaled_prior_variance = 0.2
  return(list(X=X,y=y, d=diag(t(X)%*%X), n=n,p=p,r=r,V=scaled_prior_variance * cov(y), b=beta))
}

#' @title Get lfsr per condition per CS
#' @keywords internal
mmbr_get_one_cs_lfsr = function(lfsr, alpha, sets) {
    for (i in 1:nrow(lfsr)) {
      if (i %in% sets$cs_index) {
       pos = sets$cs[[which(sets$cs_index == i)]]
       zeroed = which(!(1:nrow(lfsr) %in% pos))
       alpha[i, zeroed] = 0
     } else {
       alpha[i, ] = 0
     }
    }
    true_sign_mat = alpha * (1 - lfsr)
    pmax(1E-20, 1 - rowSums(true_sign_mat))
}

#' @title Local false sign rate (lfsr) for credible sets
#' @details This computes the lfsr of CS identified for each condition.
#' @param m a mmbr fit, the output of `mmbr::susie()`
#' @return a L by R matrix of lfsr
#' @export
mmbr_get_cs_lfsr = function(m) {
    do.call(cbind, lapply(1:dim(m$lfsr)[3], function(r) mmbr_get_one_cs_lfsr(m$lfsr[,,r], m$alpha, m$sets)))
}

#' @title Get lfsr per condition per variable
#' @importFrom stats pnorm
#' @keywords internal
mmbr_get_one_variable_lfsr = function(lfsr, alpha) {
    true_sign_mat = alpha * (1 - lfsr)
    pmax(1E-20, 1 - colSums(true_sign_mat))
}

#' @title Local false sign rate (lfsr) for variables
#' @details This computes the lfsr of variables for each condition.
#' @param m a mmbr fit, the output of `mmbr::susie()`
#' @param weighted TRUE to weight lfsr by PIP; FALSE otherwise.
#' @return a P by R matrix of lfsr
#' @export
mmbr_get_lfsr = function(m, weighted = TRUE) {
  if (weighted) alpha = m$alpha
  else alpha = matrix(1, nrow(m$alpha), ncol(m$alpha))
  do.call(cbind, lapply(1:dim(m$lfsr)[3], function(r) mmbr_get_one_variable_lfsr(m$lfsr[,,r], alpha)))
}

#' @title Make bubble heatmap to display mmbr result
#' @param m a mmbr fit, the output of `mmbr::susie()`
#' @return a plot object
#' @export
mmbr_plot = function(m, weighted_lfsr = FALSE, cs_only = TRUE, original_sumstat = FALSE) {
  if (original_sumstat) {
    if (!("bhat" %in% names(m)) || !("shat") %in% names(m))
      stop("The original summary statistics 'bhat' and 'shat' should present in input object in order to plot original summary statistics")
    if (nrow(m$bhat) != nrow(m$coef[-1,]) || nrow(m$shat) != nrow(m$coef[-1,]))
      stop(paste("Summary statistic matrix should have", nrow(m$coef[-1,]), "rows (no intercept term)."))
    bhat = m$bhat
    p = pnorm(-abs(m$bhat/m$shat))
  } else {
    bhat = m$coef[-1,]
    p = mmbr_get_lfsr(m, weighted = weighted_lfsr)
  }
  # get table of effect size estimates and PIP, for all conditions.
  table = data.frame(matrix(NA, prod(dim(p)), ncol(bhat)))
  colnames(table) = c('y', 'x', 'effect_size', 'mlog10lfsr', 'cs')
  x_names = rownames(bhat)
  y_names = colnames(bhat)
  if (is.null(x_names)) x_names = paste('variable', 1:nrow(p))
  if (is.null(y_names)) y_names = paste('condition', 1:ncol(p))
  table$y = rep(y_names, length(x_names))
  table$x = rep(x_names, each = length(y_names))
  table$effect_size = as.vector(t(bhat))
  table$mlog10lfsr = -log10(as.vector(t(p)))
  # add CS to this table.
  j = 1
  for (i in dat$sets$cs_index) {
    variables = x_names[dat$sets$cs[[j]]]
    table[which(table$x %in% variables),]$cs = i
    j = j + 1
  }
  if (cs_only) table = table[which(!is.na(table$cs)),]
  # get colors for x-axis by CS,
  xtable = unique(cbind(table$x, table$cs))
  colors = rep('black', nrow(xtable))
  for (i in unique(xtable[,2])) {
    colors[which(xtable[,2] == i)] = as.integer(i) + 2
  }
  library(ggplot2)
  p = ggplot(table) + 
    geom_point(aes(x = x, y = y, colour = effect_size , size = mlog10lfsr)) +
    scale_x_discrete(limits = unique(table$x)) + 
    scale_y_discrete(limits = unique(table$y)) + 
    scale_color_gradient2(midpoint = 0, limit = c(-max(abs(table$effect_size)), max(abs(table$effect_size))), low="#022968", mid="white", high="#800000", space="Lab") + 
    labs(size=paste0("-log10(", ifelse(original_sumstat, "p", "lfsr"), ")"), colour="Effect size") + 
    theme_minimal() + theme(text = element_text(face = "bold", size = 14), panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 15, color = colors),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())
  w = length(unique(table$x)) * 0.5
  h = length(unique(table$y)) * 0.9
  cat(paste("Suggested PDF canvas width:", w, "height:", h, "\n"))
  return(list(plot=p, width=w, height=h))
}