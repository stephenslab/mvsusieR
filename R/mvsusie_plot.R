#' @title mvSuSiE PIP and Effect Plots
#' 
#' @description Create the PIP plot and accompanying effect plot
#'   showing the effect estimates and significance of the effects for
#'   all the traits. A z-scores plot is also created when z-scores are
#'   available.
#'
#' @param fit The mvSuSiE fitted model.
#'
#' @param chr The chromosome number.
#' 
#' @param pos The positions of the genetic markers. It should have the
#'   same length as \code{fit$variable_names}.
#'
#' @param markers The names of the genetic markers (usually SNPs).
#' 
#' @param conditions The names of the conditions.
#'
#' @param poslim The range of positions to show in the PIP plot.
#'
#' @param lfsr_cutoff The significance level for lfsr. The default is
#'   0.01.
#' 
#' @param sentinel_only If \code{TRUE}, only plot the sentinel marker
#'   for each CS. If \code{FALSE}, plot all markers in each CS.
#' 
#' @param cs_plot The CSs included in the plot. The default is to show
#'   all CSs.
#' 
#' @param add_cs If \code{TRUE}, add colored dots to the top of the
#'   effect plot showing CS membership.
#' 
#' @param conditional_effect If \code{TRUE}, plot the conditional
#'   effect. If \code{FALSE}, plot the marginal effect.
#'   \code{conditional_effect = TRUE} is recommended.
#'
#' @param cs_colors The color palette for CSs.
#' 
#' @return The output includes the PIP plot, effect plot, z-scores
#'   plot (if z-scores are available in \code{fit}), and the table of
#'   effect estimates at sentinel markers.
#'
#' @examples
#' # See the "mvsusie_intro" vignette for examples.
#' 
#' @importFrom stats quantile
#' @importFrom stats median
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_size
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_line
#' @importFrom ggrepel geom_text_repel
#' @importFrom cowplot theme_cowplot
#' @importFrom cowplot plot_grid
#' 
#' @export
#'
mvsusie_plot <-
  function (fit,
            chr = 1,
            pos = seq(1,length(fit$variable_names)),
            markers = fit$variable_names,
            conditions = fit$condition_names,
            poslim = range(pos),
            lfsr_cutoff = 0.01,
            sentinel_only = TRUE,
            cs_plot = names(fit$sets$cs),
            add_cs = FALSE,
            conditional_effect = TRUE,
            cs_colors = c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00",
                "#6a3d9a", "#b15928", "#a6cee3", "#b2df8a", "#fb9a99",
                "#fdbf6f", "#cab2d6", "#ffff99", "gray", "cyan")) {

    if (!inherits(fit, "susie"))
      stop("Input argument \"fit\" should be a susie fit object, such as ",
           "the output of calling function \"mvsusie\"")
    if (length(pos) != length(fit$variable_names))
      stop("Input \"pos\" should have same length as \"fit$variable_names\"")
    if (length(markers) != length(fit$variable_names))
      stop("Input \"markers\" should have same length as ",
           "\"fit$variable_names\"")
    if (length(conditions) != length(fit$condition_names))
      stop("Input \"conditions\" should have same length as ",
           "\"fit$condition_names\"")

    # Create a data frame containing the data used for plotting.
    pdat <- data.frame("pip" = fit$pip,
                       "pos" = pos,
                       "cs"  = as.character(NA),
                       stringsAsFactors = FALSE)
    
    # Add the CS assignments to the data frame.
    #
    # WARNING: What if no identified CS?
    css <- names(fit$sets$cs) 
    for (i in css) {
      j <- fit$sets$cs[[i]]
      pdat[j, "cs"] <- i
    }

    # Create a second data frame used to plot only the points included
    # in at least one CS.
    rows    <- which(!is.na(pdat$cs))
    pdat_cs <- pdat[rows,]
    
    # Keep only the genetic markers with base-pair positions inside the
    # specified limits.
    if (!missing(poslim)) {
      rows1   <- which(pdat$pos >= poslim[1] & pdat$pos <= poslim[2])
      rows2   <- which(pdat_cs$pos >= poslim[1] & pdat_cs$pos <= poslim[2])
      pdat    <- pdat[rows1,]
      pdat_cs <- pdat_cs[rows2,]
    }
    pdat_cs$cs <- factor(pdat_cs$cs)
    css <- levels(pdat_cs$cs)

    # Reorder the CSs by position, then relabel them 1 through L.
    L          <- length(css)
    cs_pos     <- sapply(fit$sets$cs[css],function (x) median(pos[x]))
    css        <- css[order(cs_pos)]
    pdat_cs$cs <- factor(pdat_cs$cs,levels = css)
    levels(pdat_cs$cs) <- 1:L

    # Add key CS statistics to the legend (size, purity).
    cs_size <- sapply(fit$sets$cs[css],length)
    for (i in 1:L) {
      j <- css[i]
      if (cs_size[i] == 1)
        levels(pdat_cs$cs)[i] <- sprintf("%s (1 SNP)",j)
      else
        levels(pdat_cs$cs)[i] <-
          sprintf("%s (%d SNPs, %0.3f purity)",j,cs_size[j],
                  fit$sets$purity[j, "min.abs.corr"])
    }
  
    # Create a data frame containing data about the genetic markers
    # (SNPs) in the CSs (trait-specific effects, lfsr's, sentinel
    # SNPs).
    traits <- conditions
    r      <- length(traits)
    lmax   <- nrow(fit$alpha)
    fit$b1_rescaled <- fit$b1_rescaled[,-1,]
    rownames(fit$b1_rescaled) <- paste0("L",1:lmax)
    rownames(fit$single_effect_lfsr) <- paste0("L",1:lmax)
    colnames(fit$single_effect_lfsr) <- traits
    rownames(fit$alpha) <- paste0("L",1:lmax)
    effects             <- matrix(0,r,L)
    rownames(effects)   <- traits
    colnames(effects)   <- css
    effect_dat <-
      data.frame(matrix(as.numeric(NA),
                        prod(length(conditions) * length(markers)),8))
    names(effect_dat)   <- c("trait","marker","pos","effect","z","lfsr","cs",
                             "sentinel")
    effect_dat$trait    <- rep(conditions,length(markers))
    effect_dat$marker   <- rep(markers,each = length(conditions))
    effect_dat$pos      <- rep(pos,each = length(conditions))
    effect_dat$sentinel <- 0
    for (i in 1:L) {
      l <- css[i]
      j <- fit$sets$cs[[l]]
      b <- fit$b1_rescaled[l,j,]
      if (conditional_effect)
        b <- b/fit$alpha[l,j]
      marker_names <- markers[j]
      marker_idx   <- which(effect_dat$marker %in% marker_names)
      effect_dat[marker_idx,"cs"]     <- l
      effect_dat[marker_idx,"lfsr"]   <- rep(fit$single_effect_lfsr[l,],
                                              length(marker_names))
      effect_dat[marker_idx,"effect"] <- as.vector(t(b))
      if (!is.null(fit$z))
        effect_dat[marker_idx,"z"] <- as.vector(t(fit$z[j,]))
        max_idx <- which.max(fit$alpha[l,j])
      effect_dat[which(effect_dat$marker == marker_names[max_idx]),
                  "sentinel"] <- 1
      effects[,i] <- ifelse(is.null(nrow(b)),b,b[max_idx,])
    }
    effect_dat <- effect_dat[which(!is.na(effect_dat$cs)),]
    if (!missing(poslim)) {
      rows1 <- which(effect_dat$pos >= poslim[1] & effect_dat$pos <= poslim[2])
      effect_dat <- effect_dat[rows1, ]
    }
    effect_dat$marker_cs <- paste0(effect_dat$marker,"(",effect_dat$cs,")")
    pdat_sentinel <- effect_dat[which(effect_dat$sentinel == 1),]
    pdat_sentinel <- unique(pdat_sentinel[,c("marker","marker_cs","pos")])
    pdat_sentinel$pip <-
      fit$pip[match(pdat_sentinel$marker,fit$variable_names)]
    if (sentinel_only)
      effect_dat <- effect_dat[which(effect_dat$sentinel == 1),]
    if (!missing(cs_plot))
      effect_dat <- effect_dat[which(effect_dat$cs %in% cs_plot),]
    effect_dat$cs    <- factor(effect_dat$cs)
    effect_dat$trait <- factor(effect_dat$trait,traits)
    
    # Remove from the effects plot any effects that don't meet the lfsr
    # cutoff.
    rows <- which(effect_dat$lfsr < lfsr_cutoff)
    effect_dat <- effect_dat[rows,]

    # Create the PIP plot.
    pip_plot <- ggplot(pdat,aes_string(x = "pos",y = "pip")) +
      geom_point(color = "darkblue",shape = 20,size = 1.25) +
      geom_point(shape = 1,size = 1.25,stroke = 1.25,data = pdat_cs,
                 mapping = aes_string(x = "pos",y = "pip",color = "cs")) +
      geom_text_repel(data = pdat_sentinel,
                      mapping = aes_string(x="pos",y="pip",label="marker_cs"),
                      size = 2.2,segment.size = 0.35,max.overlaps = Inf,
                      min.segment.length = 0) +
      xlim(poslim[1],poslim[2]) +
      scale_color_manual(values = cs_colors) +
      guides(color = guide_legend(override.aes = list(shape = 20,
                                                      size = 1.5))) +
      labs(x = sprintf("chromosome %d position (Mb)",chr),
           y = "PIP",color = "CS") +
      theme_cowplot(font_size = 9)

    # Create the effect plot.
    if (nrow(effect_dat) > 0) {
      effect_dat$effect_sign <- factor(effect_dat$effect > 0)
      effect_dat$effect_size <- abs(effect_dat$effect)
      levels(effect_dat$effect_sign) <- c("-1","+1")
      effect_plot <- ggplot(effect_dat,
        aes_string(x = "marker",y = "trait",fill = "effect_sign",
                   size = "effect_size")) +
        geom_point(shape = 21,stroke = 0.5,color = "white") +
        scale_y_discrete(drop = FALSE) +
        scale_fill_manual(values = c("darkblue","red"),drop = FALSE) +
        scale_size(range = c(1, 5),
                   breaks = unname(quantile(effect_dat$effect_size,
                                            seq(0,1,length.out = 4)))) +
        labs(x = "",y = "",fill = "effect sign",size = "effect size") +
        guides(fill = guide_legend(override.aes = list(size = 2)),
               size = guide_legend(override.aes = list(shape=21,fill="black",
                                                       color="white",
                                                       stroke = 0.5))) +
        theme_cowplot(font_size = 9) +
        theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
              panel.grid = element_line(color = "lightgray",size = 0.3,
                                        linetype = "dotted"))

        # If requested, add colored dots to the top of the plot showing CS
        # membership.
        if (add_cs) {
          effect_dat$one <- 1
          p_cs <- ggplot(effect_dat,
                         aes_string(x = "marker",y = "one",color = "cs")) +
            geom_point(shape = 20,size = 2.5) +
            scale_x_discrete(drop = FALSE) +
            scale_color_manual(values = cs_colors) +
            labs(x = "",y = "") +
            theme_cowplot(font_size = 9) +
            theme(axis.text   = element_blank(),
                  axis.ticks  = element_blank(),
                  axis.line   = element_blank(),
                  legend.position = "top")
          effect_plot <- plot_grid(p_cs,effect_plot,nrow = 2,ncol = 1,
                                   rel_heights = c(1,3),axis = "lr",
                                   align = "v")
        }
    } else
      effect_plot <- NULL

    if (nrow(effect_dat) > 0 && all(!is.na(effect_dat$z))) {
      effect_dat$z_sign <- factor(effect_dat$z > 0)
      effect_dat$z_size <- abs(effect_dat$z)
      levels(effect_dat$z_sign) <- c("-1","+1")
      z_plot <- ggplot(effect_dat,
                       aes_string(x = "marker",y = "trait",fill = "z_sign",
                                  size = "z_size")) +
        geom_point(shape = 21,stroke = 0.5,color = "white") +
        scale_y_discrete(drop = FALSE) +
        scale_fill_manual(values = c("darkblue","red"),drop = FALSE) +
        scale_size(range = c(1,5),
                   breaks = unname(quantile(effect_dat$z_size,
                                            seq(0,1,length.out = 4)))) +
        labs(x = "",y = "",fill = "z-score sign",size = "z-score size") +
        guides(fill = guide_legend(override.aes = list(size = 2)),
               size = guide_legend(override.aes = list(shape = 21,
                                                       fill = "black",
                                                       color = "white",
                                                       stroke = 0.5))) +
        theme_cowplot(font_size = 9) +
        theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
              panel.grid = element_line(color = "lightgray",size = 0.3,
                                        linetype = "dotted"))
      if (add_cs)
        z_plot <- plot_grid(p_cs,z_plot,nrow = 2,ncol = 1,
                            rel_heights = c(1,3),axis = "lr",
                            align = "v")
    } else
      z_plot <- NULL

    # Output the (1) PIP plot, (2) effect plot, (3) z-scores plot, and
    # (4) the table of effect estimates.
    return(list(pip_plot    = pip_plot,
                effect_plot = effect_plot,
                z_plot      = z_plot,
                effects     = effects))
}
