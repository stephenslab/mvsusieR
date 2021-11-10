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
      pmax(1e-20, 1 - apply(true_sign_mat, 2, max))
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
        table[which(table$x %in% variables),]$mlog10lfsr = rep(-log10(pmax(1e-20, m$single_effect_lfsr[i,])),
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
                          low="#4575B4", mid="#FFFF0E", high="#D73027", space="Lab", na.value="white") +
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

# SuSiE model extractor
# 
#' @importFrom abind abind
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
