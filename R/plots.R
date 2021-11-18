#' @title Make bubble heatmap to display mvsusie result
#' 
#' @param m A mvsusie fit, typically the result of calling
#'   \code{\link{susie}}.
#' 
#' @return a plot object
#' 
#' @details If plot_z is TRUE, the bubble size is -log10(p value), the
#' bubble color represents effect size, the color on the x axis
#' represent CSs. If plot_z is FALSE, the bubble size is -log10(CS
#' condition specific lfsr), the bubble color represents posterior
#' effect size, the color on the x axis represent CSs, the
#' non-significant CS are removed.
#'
#' @importFrom stats pnorm
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 scale_color_gradient2
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 guide_colorbar
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 geom_rect
#' 
#' @export
#' 
mvsusie_plot = function (m, weighted_effect = FALSE, cs_only = TRUE,
                         plot_z = FALSE, pos = NULL,
                         cslfsr_threshold = 0.05) {
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
    geom_point(aes_string(x = "x",y = "y",color = "effect_size",
                          size = "mlog10lfsr")) +
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
