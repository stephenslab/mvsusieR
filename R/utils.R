#' @title chol decomposition without warning message
#' @keywords internal
muffled_chol = function(x, ...)
  withCallingHandlers(chol(x, ...),
                      warning = function(w) {
                        if (grepl("the matrix is either rank-deficient or indefinite", w$message))
                          invokeRestart("muffleWarning")
                      })

#' @title Invert a symmetric, positive definite square matrix via its Choleski decomposition
#' @keywords internal
invert_via_chol = function(x) {
  if (all(x==0)) return(x)
  return(chol2inv(muffled_chol(x)))
}

#' @title Invert SPD via triangular back-fitting
#' @keywords internal
invert_chol_tri = function(x) {
  return(t(backsolve(muffled_chol(x), diag(nrow(x)))))
}

#' @title Find trace of diag matrix
#' @keywords internal
tr = function (m) {
    if (!is.matrix(m) | (dim(m)[1] != dim(m)[2]))
        stop("Input to tr() function must be a square matrix")
    return(sum(diag(m), na.rm = TRUE))
}

#' @title Convert a list of matrices to array without losing dimension
#' @keywords internal
matlist2array = function(l) {
  if (class(l) != "list") return(l)
  l = simplify2array(l)
  if (is.null(dim(l))) {
    l = array(l, c(1,1,length(l)))
  }
  return(l)
}

#' @title compute value_j * weight_j / sum(value_j * weight_j)
#' @keywords internal
compute_softmax = function(value, weight, log = TRUE) {
    if (length(value)!=length(weight))
      stop("Values and their weights should have equal length")
    if (!log) value = log(value)
    mvalue = max(value)
    w = exp(value-mvalue)
    w_weighted = w * weight
    weighted_sum_w = sum(w_weighted)
    return(list(weights = as.vector(w_weighted / weighted_sum_w), log_sum = log(weighted_sum_w) + mvalue))
}

#' @title SuSiE model extractor
#' @importFrom abind abind
#' @keywords internal
report_susie_model = function(d, m, estimate_prior_variance = TRUE) {
    if (length(dim(m$posterior_b1[[1]])) < 2) {
      # univariate case
      b1 = t(do.call(cbind, m$posterior_b1))
      b2 = t(do.call(cbind, m$posterior_b2))
      b = colSums(b1)
    } else {
      b1 = aperm(abind::abind(m$posterior_b1,along=3), c(3,1,2))
      b2 = aperm(abind::abind(m$posterior_b2,along=3), c(3,1,2))
      if (dim(b1)[1] == 1) {
        # only one effect specified or left
        b = do.call(cbind, lapply(1:dim(b1)[3], function(i) b1[,,i]))
      } else {
        # multiple effects
        b = do.call(cbind, lapply(1:dim(b1)[3], function(i) colSums(b1[,,i])))
      }
      if (dim(b)[2] == 1) {
        b1 = b1[,,1]
        b2 = b2[,,1]
        b = as.vector(b)
      }
    }
    if (is.null(m$mixture_posterior_weights[[1]])) mixture_weights = NA
    else {
      if (length(dim(m$mixture_posterior_weights[[1]])) < 2) mixture_weights = t(do.call(cbind, m$mixture_posterior_weights))
      else mixture_weights = aperm(abind::abind(m$mixture_posterior_weights,along=3), c(3,1,2))
    }
    if (is.null(m$lfsr[[1]])) lfsr = NA
    else {
      if (length(dim(m$lfsr[[1]])) < 2) lfsr = t(do.call(cbind, m$lfsr))
      else lfsr = aperm(abind::abind(m$lfsr,along=3), c(3,1,2))
    }
    s = list(
        alpha = t(m$pip),
        b1 = b1,
        b2 = b2,
        KL = m$kl,
        lbf = m$lbf,
        V = m$prior_variance,
        sigma2 = d$residual_variance,
        elbo = m$get_objective(dump=TRUE),
        niter = m$niter,
        convergence = m$convergence,
        coef = d$rescale_coef(b),
        mixture_weights = mixture_weights,
        lfsr = lfsr
        )
    if (!is.null(m$pip_history)) s$alpha_history = m$pip_history
    if (!is.null(m$lbf_history)) s$lbf_history = m$lbf_history
    if (!is.null(m$prior_history)) s$prior_history = m$prior_history
    # FIXME: unit test for scaling issue for the fitted
    if(inherits(d, "RSSData")){
      s$fitted = d$XtX %*% b
    }else{
      s$fitted = d$compute_Xb(b)
    }
    if (is.null(dim(s$coef))) s$intercept = s$coef[1]
    else s$intercept = s$coef[1,]
    if (estimate_prior_variance) s$V = m$prior_variance
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
  condition_indicator = do.call(rbind, lapply(1:length(prior_obj$prior_variance$xUlist), function(i) as.integer(diag(prior_obj$prior_variance$xUlist[[i]]) != 0)))
  condition_pip = array(0, dim=dim(m$b1))
  for (r in 1:dim(condition_pip)[3]) {
    for (p in 1:length(condition_indicator[,r])) {
        condition_pip[,,r] = condition_pip[,,r] + m$mixture_weights[,,p] * condition_indicator[p,r]
    }
    condition_pip[,,r] = condition_pip[,,r] * m$alpha
  }
  return(condition_pip)
}

# Cannot use `unique` directly here -- for perfectly identical rows (by computation)
# due to possible numerical issues, `unique` (and `duplicated`) function reports
# that they are not identical.
almost.unique <- function(x,  tolerance = sqrt(.Machine$double.eps), ...)
{
  if (is.matrix(x)) {
    y <- round(x/tolerance, 0)
  } else {
    y <- lapply(1:length(x), function(i) round(x[[i]]/tolerance, 0))
  }
  d <- duplicated(y, ...)
  if (is.matrix(x)) x[!d,,drop=F]
  else x[!d]
}

#' @title A null progressbar, because currently `progressbar_enabled` feature does not work for `progress_bar`
#' @importFrom R6 R6Class
#' @keywords internal
null_progress_bar = R6Class('null_progress_bar', public = list(tick = function(...) {}))

#' @title check if all elements are the same in matrix of J by R, J >> R
#' @keywords internal
is_mat_common = function(mat) {
  nrow(almost.unique(mat)) == 1
}

#' @title check if all elements are the same in list
#' @keywords internal
is_list_common = function(lst) {
  length(almost.unique(lst)) == 1
}

#' @title remove duplicated columns in matrix while keeping track of what columns are removed for duplicate with what other column
#' @keywords internal
rm_collinear = function(mat, ...) {
    # `duplicated` will only work for matrix not data frame
    mat = as.matrix(mat)
    dimmat = dim(mat)
    bool_coll = almost.duplicated(mat, MARGIN = 2, ...)
    if (any(bool_coll)) {
        # these are columns to be removed
        rmvd_coll = which(bool_coll)
        # now find columns they are collinear with
        # the idea is, when using fromLast, the previously NOT duplicated columns (FALSE) will now become duplicated (TRUE)
        # then we can find these columns and use them as the columns that has some duplicated associated with them
        bool_with_coll = almost.duplicated(mat, MARGIN = 2, fromLast = TRUE, ...) & !bool_coll
        mat_with_coll = mat[, bool_with_coll, drop = FALSE]
        mat_coll = mat[, bool_coll, drop = FALSE]
        # these are columns with which the removed columns are collinear with
        # `match` will only work for data frames
        assoc_coll = which(bool_with_coll)[match(data.frame(mat_coll), data.frame(mat_with_coll))]
        rmvd_coll = cbind(assoc_coll, rmvd_coll)
        colnames(rmvd_coll) = c('associated', 'removed')
        mat = mat[, !bool_coll, drop = FALSE]
        # now generate index to recover the original
        # ...
    } else {
        rmvd_coll = NULL
    }
    attr(mat, "original_dim") = dimmat
    attr(mat, "collinear_cols") = rmvd_coll
    if (is.null(rmvd_coll)) {
        attr(mat, "collinear_counts") = NULL
    } else {
        attr(mat, "collinear_counts") = table(rmvd_coll[,'associated'])
        options(stringsAsFactors = FALSE)
        attr(mat, "collinear_counts") = cbind(as.integer(names(attr(mat, "collinear_counts"))), attr(mat, "collinear_counts") + 1)
        colnames(attr(mat, "collinear_counts")) = c('associated', 'counts')
        rownames(attr(mat, "collinear_counts")) = NULL
    }
    return(mat)
}

#' @title Reconstruct complete matrix (with duplicates) using stripped matrix and information regarding duplicate pattern in original matrix
#' @keywords internal
reconstruct_coll = function(mat, coll_cols, coll_counts, original_dim, adjust_counts=FALSE, transpose=FALSE) {
    # usage:
    # m = rm_collinear(X1)
    # X2 = reconstruct_coll(m, attr(m, 'collinear_cols'), attr(m, 'collinear_counts'), attr(m, 'original_dim'))
    # sum(X1 - X2) == 0
    get_count = function(counts, idx, adjust_counts) {
        if (!adjust_counts || ! (idx %in% counts[,"associated"])) {
            return (1)
        } else {
            print(idx)
            print(counts[,"counts"][which(counts[,"associated"] == idx)])
            return (counts[,"counts"][which(counts[,"associated"] == idx)])
        }
    }
    vectorise = FALSE
    if (is.vector(mat)) {
        vectorise = TRUE
        mat = matrix(mat, ncol = length(mat), nrow = 1)
    }
    if (transpose && !vectorise) mat = t(mat)
    # create empty matrix to the original scale
    res = matrix(NA, original_dim[1], original_dim[2])
    # first column should always be good,
    # and also duplicated columns always can be found in columns already established
    res[,1] = mat[,1] / get_count(coll_counts, 1, adjust_counts)
    i = 2
    for (j in 2:ncol(res)) {
        if (j %in% coll_cols[,'removed']) {
            # a duplicate column, just add it from before
            j0 = coll_cols[,'associated'][which(coll_cols[,'removed'] == j)]
            res[,j] = res[,j0]
        } else {
            # a new column; have to take it from the next inline in input mat
            res[,j] = mat[,i] / get_count(coll_counts, j, adjust_counts)
            i = i + 1
        }
    }
    if (transpose && !vectorise) res = t(res)
    if (vectorise) res = as.vector(res)
    return(res)
}

#' @title A simple simulation function to simulate some test data
#' @param n number of samples
#' @param p number of features
#' @param r number of conditions
#' @param s number of effect variables per condition if greater than 1; otherwise percentage of effect variables per condition
#' @param center_scale FALSE by default
#' @export
mmbr_sim1 = function(n=200,p=500,r=2,s=4,center_scale=FALSE,y_missing=NULL) {
  X = matrix(rnorm(n*p,0,1),n,p)
  if (s>=1) {
    beta = matrix(0, p, r)
    for (i in 1:r) beta[sample(1:p,s), i] = 1
  } else {
    beta = matrix(runif(p*r)>s, p, r)
  }
  y = X %*% beta + do.call(cbind, lapply(1:r, function(i) rnorm(n)))
  if (center_scale) {
    X = scale(X)
    y = t(t(y) - apply(y,2,mean))
  }
  if (!is.null(y_missing)) {
    y2 = y
    for (i in 1:nrow(y2)) {
      y2[i,runif(r) <= y_missing] = NA
    }
    y_missing = y2
  }
  scaled_prior_variance = 0.2
  return(list(X=X,y=y,y_missing=y_missing,d=diag(t(X)%*%X), n=n,p=p,r=r,V=scaled_prior_variance * cov(y),b=beta))
}

#' @title Get lfsr per condition per CS
#' @param lfsr a L by P matrix of single effect lfsr
#' @param alpha a L by P matrix of cross-condition posterior inclusion probability
#' @param sets a list of credible set output from SuSiE model
#' @keywords internal
mmbr_get_one_cs_lfsr = function(lfsr, alpha, sets) {
    # fix data dimension issue due to R's ingenuity
    if (is.null(nrow(lfsr))) lfsr = matrix(lfsr, 1, length(lfsr))
    for (i in 1:nrow(lfsr)) {
      if (i %in% sets$cs_index) {
       pos = sets$cs[[which(sets$cs_index == i)]]
       zeroed = which(!(1:ncol(lfsr) %in% pos))
       alpha[i, zeroed] = 0
       # normalize them to sum to one
       alpha[i,] = alpha[i,] / sum(alpha[i,])
     } else {
       alpha[i, ] = 0
     }
    }
    true_sign_mat = alpha * (1 - lfsr)
    pmax(0, 1 - rowSums(true_sign_mat))
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
    logp = -log10(p)
    top_snp = which(logp == max(logp, na.rm=TRUE), arr.ind = TRUE)[1]
  } else {
    bhat = m$coef[-1,]
    p = mmbr_get_lfsr(m, weighted = weighted_lfsr)
    top_snp = NULL
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
  colors = rep('black', length(table$x))
  # add CS to this table.
  if (!is.null(m$sets$cs_index)) {
    j = 1
    for (i in m$sets$cs_index) {
      variables = x_names[m$sets$cs[[j]]]
      table[which(table$x %in% variables),]$cs = i
      j = j + 1
    }
    # mark top_snp if available, as a CS by itself
    if (!is.null(top_snp)) {
        table[which(table$x == x_names[top_snp]),]$cs = 0
    }
    if (cs_only) table = table[which(!is.na(table$cs)),]
    # get colors for x-axis by CS,
    xtable = unique(cbind(table$x, table$cs))
    for (i in unique(xtable[,2])) {
      colors[which(xtable[,2] == i)] = as.integer(i) + 2
    }
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

#' @title Predict future observations or extract coefficients from susie fit
#' @param object a susie fit
#' @param newx a new value for X at which to do predictions
#' @return a matrix of predicted values for each condition
#' @details This function computes predicted values from a susie fit and a new value of X
#' @export predict.mmbr
#' @export
predict.mmbr <- function (object, newx) {
      s <- object
  for(i in 1:ncol(s$coef)){
          if(i==1){res <- s$intercept[i] + newx %*% s$coef[-1, i]} else if(i>1){
                    res <- cbind(res, s$intercept[i] + newx %*% s$coef[-1, i])
      }
    }
    return(res)
}

#' @title Compute a list of canonical covariance matrices
#' @param R an integer indicating the number of conditions
#' @param singletons a logical value indicating whether the singleton matrices are computed
#' @param hetgrid a vector of numbers between -1 and 1, each representing the off-diagonal elements of matrices with 1s on the diagonal.
#' If 0 is included, the identity matrix will be returned which corresponds to assuming effects are independent across conditions.
#' IF NULL, these matrices are not computed.
#' @return a list of canonical covariance matrices
#' @details This function computes canonical covariance matrices to be provided to mash
#' @examples
#'  mmbr:::create_cov_canonical(3)
#'  mmbr:::create_cov_canonical(3, singletons=F)
#'  mmbr:::create_cov_canonical(3, hetgrid=NULL)
#' @keywords internal
create_cov_canonical <- function(R, singletons=T, hetgrid=c(0, 0.25, 0.5, 0.75, 1)){
  mats <- list()
  s_idx <- 0
  nms <- vector()
  ###Singleton matrices
  if(singletons) {
    for(i in 1:R) {
      mats[[i]] <- matrix(0, nrow=R, ncol=R)
      mats[[i]][i, i] <- 1
      nms[i] = paste0('singleton_', i)
    }
    s_idx <- R
  }
  ###Heterogeneity matrices
  if(!is.null(hetgrid)) {
    for(j in 1:length(hetgrid)) {
      mats[[s_idx+j]] <- matrix(1, nrow=R, ncol=R)
      mats[[s_idx+j]][lower.tri(mats[[s_idx+j]], diag = FALSE)] <- hetgrid[j]
      mats[[s_idx+j]][upper.tri(mats[[s_idx+j]], diag = FALSE)] <- hetgrid[j]
      nms[s_idx+j] = paste0('shared_', j)
    }
  }
  names(mats) = nms
  return(mats)
}

#' @title Create mash prior object
#' @param fitted_g from mashr::mash
#' @param mixture_prior a list of (weights = vector(), matrices = list()) where  matrices is a list of prior matrices and have same length as weights.
#' @param sample_data a list of (X=X,Y=Y,residual_variance=residual_variance,center=T,scale=T) to allow for automatically determine canonical priors with equal weights
#' @param null_weight whether or not to add a weight for null in single effect models. By default it takes the null weight from fitted_g
#' if available. Use `null_weight = 0` to override the behavior.
#' @param use_grid expand mixture by grid values as in MASH (not necessary when prior scalar is estimated)
#' @param weights_tol filter out priors with weights smaller than weights_tol
#' @param max_mixture_len only keep the top priors by weight so that the list of mixture prior is of max_mixture_len.
#' Use `max_mixture_len=-1` to include all input weights after weights_tol filtering. Default is set to length 40.
#' @param include_indices postprocess input prior to only include conditions from this indices
#' @param ... other parameters, for mmbr:::create_cov_canonical
#' @return mash prior object for use with msusie() function
#' @details ...
#' @export
create_mash_prior = function(fitted_g = NULL, mixture_prior = NULL, sample_data = NULL,
                             null_weight = NULL, use_grid = FALSE, weights_tol = 1E-10, max_mixture_len = 40,
                             include_indices = NULL, ...) {
  if (sum(is.null(fitted_g), is.null(mixture_prior), is.null(sample_data)) != 2)
    stop("Require one and only one of fitted_g, mixture_prior and sample_data to be not NULL.")
  if (!is.null(fitted_g)) {
    # fitted_g: list(pi=pi_s, Ulist=Ulist, grid=grid, usepointmass=usepointmass)
    for (item in c('pi', 'Ulist', 'grid', 'usepointmass')) {
      if (!(item %in% names(fitted_g))) stop(paste("Cannot find", item, "in fitted_g input"))
    }
    if (fitted_g$usepointmass) {
      prior_weights = fitted_g$pi[-1]
      if (is.null(null_weight)) null_weight = fitted_g$pi[1]
    } else {
      prior_weights = mash$fitted_g$pi
    }
    return(MashInitializer$new(fitted_g$Ulist, fitted_g$grid,
                               prior_weights=prior_weights, null_weight=null_weight,
                               weights_tol=weights_tol, top_mixtures=max_mixture_len,
                               include_conditions=include_indices))
  }
  if (!is.null(mixture_prior)) {
    for (item in c('matrices')) {
      if (!(item %in% names(mixture_prior))) stop(paste("Cannot find", item, "in mixture_prior input"))
    }
    if (is.null(mixture_prior$weights)) mixture_prior$weights = rep(1/length(mixture_prior$matrices), length(mixture_prior$matrices))
    if (is.null(null_weight)) null_weight = 0
    return(MashInitializer$new(NULL, NULL, xUlist=mixture_prior$matrices, prior_weights=mixture_prior$weights,
                               null_weight=null_weight,
                               weights_tol=weights_tol, top_mixtures=max_mixture_len,
                               include_conditions=include_indices))
  }
  if (!is.null(sample_data)) {
    for (item in c('X', 'Y', 'residual_variance')) {
      if (!(item %in% names(sample_data))) stop(paste("Cannot find", item, "in sample_data input"))
    }
    if (is.null(sample_data$center)) {
      write("Assuming intercept is fitted (otherwise please set 'sample_data$center=F')", stderr())
      sample_data$center = T
    }
    if (is.null(sample_data$scale)) {
      write("Assuming X is not yet scaled and will scale X (otherwise please set 'sample_data$scale=F')", stderr())
      sample_data$scale = T
    }
    # compute canonical covariances
    Ulist = create_cov_canonical(ncol(sample_data$Y), ...)
    # compute grid
    if (use_grid) {
      d = DenseData$new(sample_data$X, sample_data$Y)
      d$standardize(sample_data$center, sample_data$scale)
      res = d$get_sumstats(diag(sample_data$residual_variance), cov2cor(sample_data$residual_variance))
      ## Use sqrt(3) giving a coarser grid than mash default in exchange for less walltime
      grid = mashr:::autoselect_grid(list(Bhat=res$bhat, Shat=res$sbhat), sqrt(3))
    } else {
      grid = 1
    }
    comp_len = length(grid) * length(Ulist)
    if (max_mixture_len<comp_len && max_mixture_len>0)
      warning(paste0('Automatically generated uniform mixture prior is of length ', comp_len, ' and is greater than currently specified max_mixture_len ', max_mixture_len, ". Please set max_mixture_len=-1 to allow using all of them (although computational speed will suffer)."))
    if (is.null(null_weight)) null_weight = 0
    return(MashInitializer$new(Ulist, grid,
                               prior_weights=NULL, null_weight=null_weight,
                               weights_tol=weights_tol, top_mixtures=max_mixture_len,
                               include_conditions=include_indices))
  }
}

#' @title Check if matrix is diag
#' @keywords internal
is_diag_mat = function(x, tol=1E-10) {
    y <- x
    diag(y) <- rep(0, nrow(y))
    return(all(abs(y) < tol))
}

#' @title Check if matrix has constant columns
#' @keywords internal
is_zero_variance <- function(x) {
  if (length(unique(x))==1) return(T)
  else return(F)
}