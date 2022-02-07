#' @title Make Bubble Heatmap to Display mvsusie Result
#' 
#' @param m A mvsusie fit, typically the result of calling
#'   \code{\link{mvsusie}}. This function only works for a mvsusie model
#'   fitted to multivariate response data (i.e., Y with more than one
#'   column) and a mash prior. mvsusie fits not satisfying these
#'   conditions will result in an error.
#'
#' @param weighted_effect If \code{weighted_effect = TRUE}, show
#'   estimated effect sizes in the original scale of X.
#'
#' @param cs_only When \code{cs_only = TRUE}, show only variants in
#'   CSs.
#'
#' @param plot_z When \code{plot_z = FALSE}, the bubble size is
#'   \eqn{-log_{10}(y)}, where \eqn{y} is the CS condition-specific
#'   lfsr, the bubble color represents posterior effect size. When
#'   \code{plot_z = TRUE}, the dots are coloured by the z-scores
#'   provided as input, and the bubble size is \eqn{-log_{10}(p)}, where
#'   \eqn{p} is the p-value. Note that \code{plot_z = TRUE} can only be
#'   used when \code{m} is an output from \code{\link{mvsusie_rss}}, or
#'   when \code{\link{mvsusie}} is called with
#'   \code{compute_univariate_zscore = TRUE}.
#'
#' @param pos A vector of indices specifying the positions of variants
#'  to show in the plot.
#'
#' @param cslfsr_threshold Effect estimates with local false sign rate
#'   (lfsr) greater than \code{cslfsr_threshold} are hidden.
#'
#' @param add_cs If \code{add_cs = TRUE}, information about the CSs is
#'   added to the top of the plot.
#' 
#' @param font_size Font size used in plot.
#' 
#' @return The return value is a list with three list elements: a
#'   ggplot object encoding the generated plot, and \dQuote{height} and
#'   \dQuote{width} giving the suggested plot dimensions.
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
#' @importFrom ggplot2 scale_color_manual
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
#' @importFrom cowplot plot_grid
#'
#' @export
#' 
mvsusie_plot = function (m, weighted_effect = FALSE, cs_only = TRUE,
                         plot_z = FALSE, pos = NULL, cslfsr_threshold = 0.05,
                         add_cs = TRUE, font_size = 12) {

  cs_colors = c(
    "#FF7F00", # orange
    "skyblue2",
    "green1",
    "#6A3D9A", # purple
    "#FB9A99", # lt pink
    "dodgerblue2",
    "green4",
    "gold1",
    "palegreen2",
    "#CAB2D6", # lt purple
    "#FDBF6F", # lt orange
    "gray70", "khaki2",
    "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
    "darkturquoise", "yellow4", "yellow3",
    "darkorange4", "brown"
  )
  
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
  } else {
    if (weighted_effect)
      effects = m$coef[-1,]
    else
      effects = colSums(m$b1)
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
  table$one         = 1
  table$mlog10lfsr  = NA
  table$effect_size = NA
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
        y = pmax(1e-20,m$single_effect_lfsr[i,])
        y[y > cslfsr_threshold] = NA
        table[which(table$x %in% variables),"mlog10lfsr"] =
          rep(-log10(y),length(variables))
        idx = which((table$x %in% variables) & (table$y %in% condition_sig))
        table[idx,"effect_size"] = effects[idx]
      }
      j = j + 1
    }
    if (cs_only)
      table = table[which(!is.na(table$cs)),]
  }
  rowidx            = which(table$x %in% x_names[pos])
  table             = table[rowidx,]
  a                 = min(table$effect_size,na.rm = TRUE) - 1e-8
  b                 = max(table$effect_size,na.rm = TRUE)
  table$x           = factor(table$x)
  xlabels           = levels(table$x)
  xlabels[tapply(table$cs,table$x,
                 function (x) all(is.na(x)))] <- ""
  table$effect_size = cut(table$effect_size,
                          breaks = c(seq(a,-1e-8,length.out = 4),
                                   c(seq(1e-8,b,length.out = 4))))
  table$cs          = factor(table$cs)
  
  p = ggplot(table) +
    geom_point(mapping = aes_string(x = "x",y = "y",fill = "effect_size",
                                    size = "mlog10lfsr"),
               shape = 21,color = "white",na.rm = TRUE) +
    scale_x_discrete(limits = unique(table$x),labels = xlabels,drop = FALSE) +
    scale_y_discrete(limits = unique(table$y),drop = FALSE) + 
    scale_radius(range = c(1,10)) +
        
    # Colors obtained from colorbrewer2.org.
    scale_fill_manual(values = c("darkblue","#0571b0","#92c5de","gainsboro",
                                 "#f4a582","#ca0020","firebrick"),
                      na.value = "white",drop = FALSE) +
    labs(size = paste0("-log10(",ifelse(plot_z,"p","CS lfsr"),")"),
         fill = ifelse(plot_z, 'z scores', 'NCP')) +
    guides(
      size = guide_legend(order = 1, override.aes = list(color = "black",fill = "black")),
      fill = guide_legend(order = 2, override.aes = list(size = 3))) +
    theme_cowplot(font_size = font_size) +
    theme(panel.grid   = element_blank(),
          axis.text.x  = element_text(angle = 45,vjust = 1,hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())

  # If requested, add colored dots to the top of the plot showing CS
  # membership.
  if (add_cs) {
    rows <- which(!is.na(table$cs))
    p_cs = ggplot(table[rows,],aes_string(x = "x",y = "one",color = "cs")) +
      geom_point(shape = 20,size = 2.5) +
      scale_x_discrete(drop = FALSE) +
      scale_color_manual(values = cs_colors) +
      labs(x = "",y = "") +
      theme_cowplot(font_size = font_size) +
      theme(axis.text   = element_blank(),
            axis.ticks  = element_blank(),
            axis.line   = element_blank())
    p = plot_grid(p_cs,p,nrow = 2,ncol = 1,rel_heights = c(1,3),
                  axis = "lr",align = "v")
  }
    
  w = 3 + 0.6 * length(unique(table$x))
  h = 0.7 * length(unique(table$y))
  cat("Suggested PDF canvas width:",w,"height:",h,"\n")
  return(list(plot = p,width = w,height = h))
}
