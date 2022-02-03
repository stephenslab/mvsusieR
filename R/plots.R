#' @title Make Bubble Heatmap to Display mvsusie Result
#' 
#' @param m A mvsusie fit, typically the result of calling
#'   \code{\link{mvsusie}}. This function only works for a mvsusie model
#'   fitted to multivariate response data (i.e., Y with more than one
#'   column) and a mash prior. mvsusie fits not satisfying these
#'   conditions will result in an error.
#'
#' @param weighted_effect Describe input argument "weighted_effect"
#'   here.
#'
#' @param cs_only Describe input argument "cs_only" here.
#'
#' @param plot_z If \code{plot_z = TRUE}, the bubble size is
#'   \eqn{-log_{10}(p)}, where \eqn{p} is the p-valule. The bubble color
#'   represents effect size, the color on the x-axis represent CSs. If
#'   \code{plot_z = FALSE}, the bubble size is \eqn{-log_{10}(y)}, where
#'   \eqn{y} is the CS condition-specific lfsr, the bubble color
#'   represents posterior effect size, and the color on the x-axis
#'   represent CSs, with the non-significant CS removed.
#'
#' @param pos Describe input argument "pos" here.
#'
#' @param cslfsr_threshold Describe input argument "cslfsr_threshold"
#'   here.
#'
#' @param font_size Font size used in plot.
#' 
#' @return Describe output here.
#'
#' @examples
#' n = 500
#' p = 1000
#' true_eff = 2
#' X = matrix(sample(0:2,size = n*p,replace = TRUE),nrow = n,ncol = p)
#' beta1 = rep(0,p)
#' beta2 = rep(0,p)
#' beta3 = rep(0,p)
#' beta1[1:true_eff] = runif(true_eff)
#' beta2[1:true_eff] = runif(true_eff)
#' beta3[1:true_eff] = runif(true_eff)
#' y1 = X %*% beta1 + rnorm(n)
#' y2 = X %*% beta2 + rnorm(n)
#' y3 = X %*% beta3 + rnorm(n)
#' Y = cbind(y1,y2,y3)
#' prior =
#'   create_mash_prior(sample_data = list(X=X,Y=Y,residual_variance=cov(Y)),
#'                     max_mixture_len = -1)
#' fit = mvsusie(X,Y,prior_variance = prior,compute_univariate_zscore = TRUE)
#' mvsusie_plot(fit,plot_z = FALSE)$plot
#' mvsusie_plot(fit,plot_z = TRUE)$plot
#'
#' @importFrom stats pnorm
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_radius
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 guide_colorbar
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 geom_rect
#' @importFrom cowplot theme_cowplot
#'
#' @export
#' 
mvsusie_plot = function (m, weighted_effect = FALSE, cs_only = TRUE,
                         plot_z = FALSE, pos = NULL, cslfsr_threshold = 0.05,
                         font_size = 12) {

  # The susie fit should be for multivariate Y with mash prior.
  if (!inherits(m,"susie"))
    stop("Input argument \"m\" should be a susie fit object, such as the ",
         "output of calling function \"mvsusie\"")
  
  if (plot_z) {
    if (!is.element("z",names(m)))
      stop("Cannot find the z-score summary statistics")
    if (nrow(m$z) + 1 != nrow(m$coef))
      stop(paste("z-score matrix should have",nrow(m$z) + 1,"rows"))
    effects = m$z
    p       = 2*pnorm(-abs(m$z))
    logp    = -log10(p)
    top_snp = which(logp == max(logp,na.rm = TRUE),arr.ind = TRUE)[1]
  } else {
    if (weighted_effect)
      effects = m$coef[-1,]
    else
      effects = colSums(m$b1)
    top_snp = NULL
  }
  if (is.null(pos))
    pos = 1:nrow(effects)
  else if (!all(pos %in% seq(1,nrow(effects))))
    stop("Provided position is outside the range of variables")
  x_names = m$variable_names
  y_names = m$condition_names
  if (is.null(x_names))
    x_names = paste("variable",1:nrow(effects))
  if (is.null(y_names))
    y_names = paste("condition",1:ncol(effects))
  table             = data.frame(matrix(as.numeric(NA),prod(dim(effects)),5))
  colnames(table)   = c("y","x","effect_size","mlog10lfsr","cs")
  table$y           = rep(y_names,length(x_names))
  table$x           = rep(x_names,each = length(y_names))
  table$mlog10lfsr  = 0
  table$effect_size = 0
  if (plot_z) {
    table$mlog10lfsr  = as.vector(t(logp))
    table$effect_size = as.vector(t(effects))
  }
  
  # Add CS to this table.
  if (!is.null(m$sets$cs_index)) {
    if (!plot_z)
      effects = as.vector(t(effects))
    j = 1
    for (i in m$sets$cs_index) {
      condition_idx = which(m$single_effect_lfsr[i,] < cslfsr_threshold)
      condition_sig = y_names[condition_idx]
      variables = x_names[m$sets$cs[[j]]]
      table[which(table$x %in% variables),]$cs = i
      if (!plot_z) {
        table[which(table$x %in% variables),]$mlog10lfsr =
          rep(-log10(pmax(1e-20,m$single_effect_lfsr[i,])),length(variables))
        idx = which((table$x %in% variables) & (table$y %in% condition_sig))
        table[idx,]$effect_size = effects[idx]
      }
      j = j + 1
    }
    if (cs_only)
      table = table[which(!is.na(table$cs)),]
  }
  
  rowidx            = which(table$x %in% x_names[pos])
  table             = table[rowidx,]
  a                 = min(table$effect_size)
  b                 = max(table$effect_size)
  table$effect_size = cut(table$effect_size,
                          breaks = c(seq(a,-1e-8,length.out = 4),
                                   c(seq(1e-8,b,length.out = 4))))
  table$x           = factor(table$x)
  xlabels           = levels(table$x)
  
  # cs_colors = unique(cbind(table$x,table$cs,table$color))[,3]
  p = ggplot(table) +
    geom_point(mapping = aes_string(x = "x",y = "y",fill = "effect_size",
                                    size = "mlog10lfsr"),
               shape = 21,color = "white") +
    scale_x_discrete(limits = unique(table$x),labels = xlabels) +
    scale_y_discrete(limits = unique(table$y)) + 
    scale_radius(range = c(1,10)) +
    # Colors obtained from colorbrewer2.org.
    scale_fill_manual(values = c("darkblue","#0571b0","#92c5de","gainsboro",
                                 "#f4a582","#ca0020","firebrick"),
                      na.value = "gainsboro",drop = FALSE) +
    labs(size = paste0("-log10(",ifelse(plot_z,"p","CS lfsr"),")")) +
    guides(size = guide_legend(override.aes = list(color = "black",
                                                   fill = "black"))) +
    theme_cowplot(font_size = font_size) +
    theme(panel.grid   = element_blank(),
          axis.text.x  = element_text(angle = 45,vjust = 1,hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  if (!is.null(top_snp))
    if (top_snp %in% unlist(m$sets$cs)) {
      xnode = which(unique(table$x) == x_names[top_snp])
      p = p + geom_rect(data = data.frame(ymin = 0.5,
                                          ymax = length(unique(table$y)) + 0.5,
                                          xmin = xnode - 0.5,
                                          xmax = xnode + 0.5),
                        mapping = aes_string(ymin = "ymin",
                                             ymax = "ymax",
                                             xmin = "xmin",
                                             xmax = "xmax"),
                        color = "black",fill = "white",alpha = 0)
    }
  w = 3 + 0.6 * length(unique(table$x))
  h = 0.7 * length(unique(table$y))
  cat("Suggested PDF canvas width:",w,"height:",h,"\n")
  return(list(plot = p,width = w,height = h))
}
