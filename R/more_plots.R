#' @title mvSuSiE PIP and Effect Plots
#' 
#' @description Create the PIP plot and accompanying effect plot
#'  showing the effect estimates and significance of the effects for
#'  all the traits.
#'
#' @details The colors from colorbrewer2.org.
#'
#' @param fit Describe the fit input here.
#'
#' @param pos Describe the pos input here.
#'
#' @param markers Describe the markers input here.
#'
#' @param chr Describe the chr input here.
#'
#' @param poslim Describe the poslim input here.
#'
#' @param conditions Describe the conditions input here.
#'
#' @param lfsr_cutoff Describe the lfsr_cutoff input here.
#'
#' @param lfsr_breaks Describe the lfsr_breaks input here.
#'
#' @param cs_colors Describe the cs_colors input here.
#' 
#' @return Describe the return value here.
#'
#' @importFrom stats quantile
#' @importFrom stats median
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_alpha_manual
#' @importFrom ggplot2 scale_size
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_line
#' @importFrom ggrepel geom_text_repel
#' @importFrom cowplot theme_cowplot
#' 
#' @export
#'
mvsusie_plot_simple <-
  function (fit,
            pos = seq(1,length(fit$variable_names)),
            markers = fit$variable_names,
            chr = 1,
            poslim = range(pos),
            conditions,
            lfsr_cutoff = 0.01,
            lfsr_breaks = c(-Inf,1e-15,1e-8,1e-4,0.05,Inf),
            cs_colors = c("#1f78b4","#33a02c","#e31a1c","#ff7f00",
                          "#6a3d9a","#b15928","#a6cee3","#b2df8a",
                          "#fb9a99","#fdbf6f","#cab2d6","#ffff99",
                          "gray","cyan")) {

  # Create a data frame containing the data used for plotting.
  pdat <- data.frame("pip" = fit$pip,
                     "pos" = pos,
                     "cs"  = as.character(NA),
                     stringsAsFactors = FALSE)

  # Add the CS assignments to the data frame.
  css <- names(fit$sets$cs)
  for (i in css) {
    j <- fit$sets$cs[[i]]
    pdat[j,"cs"] <- i
  }

  # Create a second data frame used to plot only the points included
  # in at least one CS.
  rows <- which(!is.na(pdat$cs))
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
  css        <- levels(pdat_cs$cs)
  
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
      levels(pdat_cs$cs)[i] <- sprintf("%d (1 SNP)",i)
    else
      levels(pdat_cs$cs)[i] <-
        sprintf("%d (%d SNPs, %0.3f purity)",i,cs_size[j],
                fit$sets$purity[j,"min.abs.corr"])
  }
  
  # Create two more data frames containing data about the "sentinel"
  # SNPs only: pdat_sentinel contains data about the positions of the
  # sentinel SNPs; pdat_effets contains data about the trait-wise
  # effects and lfsrs.
  traits <- fit$condition_names
  n      <- length(traits)
  lmax   <- nrow(fit$alpha)
  pdat_sentinel <- data.frame("pos"    = rep(0,L),
                              "pip"    = rep(0,L),
                              "lfsr"   = rep(0,L),
                              "marker" = rep("",L),
                              stringsAsFactors = FALSE)
  pdat_effects <- data.frame("trait"     = rep(traits,times = L),
                             "cs"        = rep(1:L,each = n),
                             "coef_sign" = 0,
                             "coef_size" = 0,
                             "lfsr"      = 0,
                             stringsAsFactors = FALSE)
  fit$b1_rescaled <- fit$b1_rescaled[,-1,]
  effects <- matrix(0,n,L)
  rownames(effects) <- traits
  colnames(effects) <- 1:L
  rownames(fit$b1_rescaled)        <- paste0("L",1:lmax)
  rownames(fit$single_effect_lfsr) <- paste0("L",1:lmax)
  colnames(fit$single_effect_lfsr) <- traits
  for (i in 1:L) {
    l    <- css[i]
    j    <- fit$sets$cs[[l]]
    j    <- j[which.max(fit$pip[j])]
    rows <- which(pdat_effects$cs == i)
    b    <- fit$b1_rescaled[l,j,]
    pdat_sentinel[i,"pos"]         <- pos[j]
    pdat_sentinel[i,"pip"]         <- fit$pip[j]
    pdat_sentinel[i,"marker"]      <- sprintf("%s (%d)",markers[j],i)
    pdat_effects[rows,"coef_size"] <- abs(b)
    pdat_effects[rows,"coef_sign"] <- (b > 0)
    pdat_effects[rows,"lfsr"]      <- fit$single_effect_lfsr[l,]
    effects[,i]                    <- b
  }
  if (!missing(conditions))
      traits <- conditions
  pdat_effects$cs         <- factor(pdat_effects$cs)
  pdat_effects$trait      <- factor(pdat_effects$trait,traits)
  pdat_effects$coef_sign  <- factor(pdat_effects$coef_sign)
  levels(pdat_effects$cs) <- pdat_sentinel$marker
  levels(pdat_effects$coef_sign) <- c("-1","+1")

  # Remove from the effects plot any effects that don't meet the lfsr
  # cutoff.
  rows <- which(pdat_effects$lfsr < lfsr_cutoff)
  pdat_effects <- pdat_effects[rows,]
  pdat_effects$lfsr <- cut(pdat_effects$lfsr,lfsr_breaks)
  
  # Create the PIP plot.
  pip_plot <- ggplot(pdat,aes_string(x = "pos",y = "pip")) +
    geom_point(color = "darkblue",shape = 20,size = 1.25) +
    geom_point(shape = 1,size = 1.25,stroke = 1.25,data = pdat_cs,
               mapping = aes_string(x = "pos",y = "pip",color = "cs")) +
    geom_text_repel(data = pdat_sentinel,
                    mapping = aes_string(x = "pos",y = "pip",label = "marker"),
                    size = 2.2,segment.size = 0.35,max.overlaps = Inf,
                    min.segment.length = 0) +
    xlim(poslim[1],poslim[2]) +
    scale_color_manual(values = cs_colors) +
    guides(color = guide_legend(override.aes = list(shape = 20,size = 1.5))) +
    labs(x = sprintf("chromosome %d position (Mb)",chr),
         y = "PIP",color = "CS") +
    theme_cowplot(font_size = 9)

  # Create the effect plot. 
  effect_plot <- ggplot(pdat_effects,
                        aes_string(x = "cs",y = "trait",fill = "coef_sign",
                                   size = "coef_size",alpha = "lfsr")) +
    geom_point(shape = 21,stroke = 0.5,color = "white") +
    scale_y_discrete(drop = FALSE) + 
    scale_fill_manual(values = c("darkblue","red"),drop = FALSE) +
    scale_alpha_manual(values = c(0.95,0.8,0.65,0.5,0.05),drop = FALSE) +
    scale_size(range = c(1,5),
               breaks = unname(quantile(pdat_effects$coef_size,
                                        seq(0,1,length.out = 4)))) +
    labs(x = "",y = "",fill = "sign",size = "size") +
    guides(alpha = guide_legend(override.aes = list(shape = 21,color = "white",
                                                    fill = "black",size = 2)),
           fill = guide_legend(override.aes = list(size = 2,alpha = 0.95)),
           size = guide_legend(override.aes = list(shape = 21,fill = "black",
                                                   color = "white",
                                                   stroke = 0.5,
                                                   alpha = 0.95))) +
    theme_cowplot(font_size = 9) +
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
          panel.grid = element_line(color = "lightgray",linetype = "dotted",
                                    size = 0.3))

  # Output the (1) PIP plot, (2) effect plot and (3) the table of
  # effect estimates.
  return(list(pip_plot     = pip_plot,
              effect_plot  = effect_plot,
              effects      = effects))
}
