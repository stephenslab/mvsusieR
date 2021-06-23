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
  if (all(x==0)) return(list(inv=x, rank=0))
  return(list(inv=chol2inv(muffled_chol(x)), rank=nrow(x)))
}

#' @title Invert SPD via triangular back-fitting
#' @keywords internal
invert_chol_tri = function(x) {
  return(list(inv=t(backsolve(muffled_chol(x), diag(nrow(x)))), rank=nrow(x)))
}

#' @title Pseudo inverse of matrix
#' @keywords internal
pseudo_inverse = function(x, tol=sqrt(.Machine$double.eps)){
  xsvd <- svd(x)
  Positive <- xsvd$d > max(tol * xsvd$d[1L], 0)
  if(all(Positive)){
    xinv <- xsvd$v %*% (1/xsvd$d * t(xsvd$u))
  }else{
    xinv <- xsvd$v[, Positive, drop = FALSE] %*% ((1/xsvd$d[Positive]) *
                                            t(xsvd$u[, Positive, drop = FALSE]))
  }
  return(list(inv = xinv, rank = sum(Positive)))
}

#' @title Check if x is diagonal matrix
#' @keywords internal
isDiagonal = function(x, tol=sqrt(.Machine$double.eps)){
  if(is.matrix(x)){
    diag(x) <- rep(0, nrow(x))
    return(all(abs(x) < tol))
  }else{
    return(TRUE)
  }
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
    if (is.null(m$clfsr[[1]])) clfsr = NA
    else {
      if (length(dim(m$clfsr[[1]])) < 2) clfsr = t(do.call(cbind, m$clfsr))
      else clfsr = aperm(abind::abind(m$clfsr,along=3), c(3,1,2))
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
        conditional_lfsr = clfsr,
        lfsr = mvsusie_get_lfsr(clfsr, t(m$pip)),
        single_effect_lfsr = mvsusie_single_effect_lfsr(clfsr, t(m$pip))
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

#' @title Compute condition specific posterior inclusion probability.
#' @description This is only relevant when canonical priors are used
#' @param m M&M model
#' @param prior_obj prior mixture object
#' @return P by R matrix of PIP per condition
#' @keywords internal
mvsusie_get_pip_per_condition = function(m, prior_obj) {
  condition_pip = mvsusie_get_alpha_per_condition(m, prior_obj)
  return(do.call(cbind, lapply(1:dim(condition_pip)[3], function(r) apply(condition_pip[,,r], 2, function(x) 1-prod(1-x)))))
}

#' @title Compute condition specific posterior inclusion probability per effect
#' @keywords internal
mvsusie_get_alpha_per_condition = function(m, prior_obj) {
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

#' @title `duplicated` function with a tolerance
#' @keywords internal
almost.duplicated <- function(x, tolerance = sqrt(.Machine$double.eps), ...)
{
  y <- round(x/tolerance, 0)
  duplicated(y, ...)
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
mvsusie_sim1 = function(n=200,p=500,r=2,s=4,center_scale=FALSE,y_missing=NULL) {
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

#' @title Local false sign rate (lfsr) for single effects
#' @details This computes the lfsr of single effects for each condition.
#' @param alpha L by P matrix
#' @param clfsr L by P by R conditonal lfsr
#' @return a L by R matrix of lfsr
#' @export
mvsusie_single_effect_lfsr = function(clfsr, alpha) {
  if(!is.array(clfsr) && is.na(clfsr)){
    return(NA)
  }else{
    return(do.call(cbind, lapply(1:dim(clfsr)[3], function(r){
      clfsrr = clfsr[,,r]
      if (is.null(nrow(clfsrr))) clfsrr = matrix(clfsrr, 1, length(clfsrr))
      pmax(0, rowSums(alpha * clfsrr))
    })))
  }
}

#' @title Local false sign rate (lfsr) for variables
#' @details This computes the lfsr of variables for each condition.
#' @param alpha L by P matrix
#' @param clfsr L by P by R conditonal lfsr
#' @param weighted TRUE to weight lfsr by PIP; FALSE otherwise.
#' @return a P by R matrix of lfsr
#' @export
mvsusie_get_lfsr = function(clfsr, alpha, weighted = TRUE) {
  if(!is.array(clfsr) && is.na(clfsr)){
    return(NA)
  }else{
    if (weighted) alpha = alpha
    else alpha = matrix(1, nrow(alpha), ncol(alpha))
    return(do.call(cbind, lapply(1:dim(clfsr)[3], function(r){
      true_sign_mat = alpha * (1 - clfsr[,,r])
      pmax(1E-20, 1 - apply(true_sign_mat, 2, max))
    })))
  }
}

#' @title Make bubble heatmap to display mvsusie result
#' @param m a mvsusie fit, the output of `mvsusieR::susie()`
#' @return a plot object
#' @details If plot_z is TRUE, the bubble size is -log10(p value), 
#'   the bubble color represents effect size, the color on the x axis represent
#'   CSs. If plot_z is FALSE, the bubble size is -log10(CS condition specific lfsr), 
#'   the bubble color represents posterior effect size, the color on the x axis represent
#'   CSs, the non-significant CS are removed. 
#' @export
mvsusie_plot = function(m, weighted_effect = FALSE, cs_only = TRUE,
                        plot_z = FALSE, pos = NULL, cslfsr_threshold = 0.05) {
  colors = c(
    "dodgerblue2",
    "green4",
    "#6A3D9A", # purple
    "#FF7F00", # orange
    "gold1",
    "skyblue2", "#FB9A99", # lt pink
    "palegreen2",
    "#CAB2D6", # lt purple
    "#FDBF6F", # lt orange
    "gray70", "khaki2",
    "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
    "darkturquoise", "green1", "yellow4", "yellow3",
    "darkorange4", "brown"
  )
  if (plot_z) {
    if (!("z" %in% names(m)))
      stop("Cannot find the z score summary statistics.")
    if (nrow(m$z) != nrow(m$coef[-1,]))
      stop(paste("z score matrix should have", nrow(m$coef[-1,]), "rows (no intercept term)."))
    effects = m$z
    p = pnorm(-abs(m$z))*2
    logp = -log10(p)
    top_snp = which(logp == max(logp, na.rm=TRUE), arr.ind = TRUE)[1]
  } else {
    if (weighted_effect) effects = m$coef[-1,]
    else effects = colSums(m$b1, dims=1)
    # if (all(is.na(m$lfsr))) stop("Cannot make bubble plot without lfsr information (currently only implemented for mixture prior)")
    # p = m$lfsr
    top_snp = NULL
  }
  if(is.null(pos)){
    pos = 1:nrow(effects)
  }else{
    if (!all(pos %in% 1:nrow(effects))) 
      stop("Provided position is outside the range of variables")
  }
  
  x_names = m$variable_names
  y_names = m$condition_names
  if (is.null(x_names)) x_names = paste('variable', 1:nrow(effects))
  if (is.null(y_names)) y_names = paste('condition', 1:ncol(effects))
  
  table = data.frame(matrix(NA, prod(dim(effects)), 6))
  colnames(table) = c('y', 'x', 'effect_size', 'mlog10lfsr', 'cs', 'color')
  table$y = rep(y_names, length(x_names))
  table$x = rep(x_names, each = length(y_names))
  table$color = 'black'
  if(plot_z){
    table$mlog10lfsr = as.vector(t(logp))
    table$effect_size = as.vector(t(effects))
  }
  
  # add CS to this table.
  if (!is.null(m$sets$cs_index)) {
    if(plot_z == FALSE){
      effects = as.vector(t(effects))
    }
    j = 1
    for (i in m$sets$cs_index) {
      condition_idx = which(m$single_effect_lfsr[i,] < cslfsr_threshold)
      condition_sig = y_names[condition_idx]
      variables = x_names[m$sets$cs[[j]]]
      table[which(table$x %in% variables),]$cs = i
      table[which(table$x %in% variables),]$color = colors[i %% length(colors)]
      if(plot_z == FALSE){
        table[which(table$x %in% variables),]$mlog10lfsr = rep(-log10(pmax(1E-20, m$single_effect_lfsr[i,])),
                                                               length(variables))
        idx = which((table$x %in% variables) & (table$y %in% condition_sig))
        table[idx,]$effect_size = effects[idx]
      }
      j = j + 1
    }
    if (cs_only){
      table = table[which(!is.na(table$cs)),]
    }
  }
  
  rowidx = which(table$x %in% x_names[pos])
  table = table[rowidx,]
  
  cs_colors = unique(cbind(table$x, table$cs, table$color))[,3]
  
  p = ggplot(table) +
    geom_point(aes(x = x, y = y, colour = effect_size, size = mlog10lfsr)) +
    scale_x_discrete(limits = unique(table$x)) +
    scale_y_discrete(limits = unique(table$y)) + 
    scale_color_gradient2(midpoint = 0, limit = c(-max(abs(table$effect_size), na.rm=TRUE), 
                                                  max(abs(table$effect_size), na.rm=TRUE)),
                          low="#4575B4", mid="#FFFFBF", high="#D73027", space="Lab", na.value="white") +
    labs(size=paste0("-log10(", ifelse(plot_z, "p", "CS lfsr"), ")"), colour=ifelse(plot_z, "z-score", "Effect size")) +
    guides(size = guide_legend(order = 1), colour = guide_colorbar(order = 2)) +
    theme_minimal() + theme(text = element_text(face = "bold", size = 14), panel.grid = element_blank(),
                            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 15, color = cs_colors),
                            axis.text.y = element_text(size = 15, color = "black"),
                            axis.title.x = element_blank(),
                            axis.title.y = element_blank())
  
  if (!is.null(top_snp)) {
    if(top_snp %in% unlist(m$sets$cs)){
      xnode = which(unique(table$x) == x_names[top_snp])
      
      p = p + geom_rect(aes(ymin=1-0.5,
                            ymax=length(unique(table$y)) + 0.5,
                            xmin=xnode-0.5,
                            xmax=xnode+0.5), color="black", alpha=0, fill = 'white')
    }
  }
  w = length(unique(table$x)) * 0.6 + 3
  h = length(unique(table$y)) * 0.7
  cat(paste("Suggested PDF canvas width:", w, "height:", h, "\n"))
  return(list(plot=p, width=w, height=h))
}


#' @title Predict future observations or extract coefficients from susie fit
#' @param object a susie fit
#' @param newx a new value for X at which to do predictions
#' @return a matrix of predicted values for each condition
#' @details This function computes predicted values from a susie fit and a new value of X
#' @export predict.mvsusie
#' @export
predict.mvsusie <- function (object, newx) {
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
#'  mvsusieR:::create_cov_canonical(3)
#'  mvsusieR:::create_cov_canonical(3, singletons=F)
#'  mvsusieR:::create_cov_canonical(3, hetgrid=NULL)
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
#' Use `max_mixture_len=-1` to include all input weights after weights_tol filtering. Default is set to use all input prior matrices.
#' @param include_indices postprocess input prior to only include conditions from this indices
#' @param ... other parameters, for mvsusieR:::create_cov_canonical
#' @return mash prior object for use with mvsusie() function
#' @details ...
#' @export
create_mash_prior = function(fitted_g = NULL, mixture_prior = NULL, sample_data = NULL,
                             null_weight = NULL, use_grid = FALSE, weights_tol = 1E-10, max_mixture_len = -1,
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
      #message("Assuming intercept is fitted (otherwise please set 'sample_data$center=F')", stderr())
      sample_data$center = T
    }
    if (is.null(sample_data$scale)) {
      #message("Assuming X is not yet scaled and will scale X (otherwise please set 'sample_data$scale=F')", stderr())
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

#'@title Scale prior matrix
#' @keywords internal
scale_covariance <- function(mat, sigma) {
  # faster way to implement diag(sigma) %*% mat %*% diag(sigma)
  t(mat * sigma) * sigma
}

#' @title Check if input is numeric matrix
#' @keywords internal
is_numeric_matrix <- function(X, name) {
  if (!((is.double(X) || is.integer(X)) & is.matrix(X)))
    stop(paste("Input", name, "must be a numeric matrix."))
  if (any(is.na(X))) {
    stop(paste("Input", name, "must not contain any missing values."))
  }
  if (any(dim(X) == 0))
    stop(paste("Input", name, "dimension is invalid."))
}
