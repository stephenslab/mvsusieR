#' @title mvSuSiE PIP and Effect Plots
#' 
#' @description Create the PIP plot and accompanying effect plot
#'  showing the effect estimates and significance of the effects for
#'  all the traits. Optionally, a z-scores plot is also created.
#'
#' @param fit The mvSuSiE fitted model.
#'
#' @param chr The chromosome number.
#' 
#' @param pos The variant position. It should have the same length as
#'   \code{fit$variable_names}.
#'
#' @param markers The variant names.
#' 
#' @param conditions The condition names.
#'
#' @param poslim The range of position in the plot.
#'
#' @param lfsr_cutoff The significant level for lfsr. The default is 0.01.
#' 
#' @param sentinel_only If TRUE, only plot the sentinel variant for each CS.
#' 
#' @param cs_plot The CSs included in the plot. The default is all CSs.
#' 
#' @param add_cs If TRUE, add colored dots to the top of the plot showing
#'   CS membership.
#' 
#' @param conditional_effect If TRUE, plot the conditional effect.
#'
#' @param cs_colors The color palette for CSs.
#' 
#' @return The output includes the PIP plot, effect plot, z-scores
#'   plot (if z scores are available in \code{fit}), and the table of
#'   effect estimates at sentinel variants.
#'
#' @examples
#' # See the "mvsusie_intro" vignette for examples.
#' 
#' @importFrom stats quantile
#' @importFrom stats median
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
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
  function (fit, chr = 1, pos = seq(1, length(fit$variable_names)),
            markers = fit$variable_names, conditions = fit$condition_names,
            poslim = range(pos), lfsr_cutoff = 0.01, sentinel_only = TRUE,
            cs_plot = names(fit$sets$cs), add_cs = FALSE,
            conditional_effect = TRUE,
            cs_colors = c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00",
                "#6a3d9a", "#b15928", "#a6cee3", "#b2df8a", "#fb9a99",
                "#fdbf6f", "#cab2d6", "#ffff99", "gray", "cyan")){

    if (!inherits(fit, "susie"))
      stop("Input argument \"fit\" should be a susie fit object, such as the ",
           "output of calling function \"mvsusie\"")

    if (length(pos) != length(fit$variable_names))
      stop("Input \"pos\" should have same length as \"fit$variable_names\"")
    if (length(markers) != length(fit$variable_names))
      stop("Input \"markers\" should have same length as \"fit$variable_names\"")
    if (length(conditions) != length(fit$condition_names))
      stop("Input \"conditions\" should have same length as \"fit$condition_names\"")

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
  
    # Create a data frame containing data about the variants
    # in the CSs (trait-specific effects, lfsrs, sentinel variant)
    traits <- conditions
    r      <- length(traits)
    lmax   <- nrow(fit$alpha)
    fit$b1_rescaled <- fit$b1_rescaled[, -1, ]
    rownames(fit$b1_rescaled) <- paste0("L", 1:lmax)
    rownames(fit$single_effect_lfsr) <- paste0("L", 1:lmax)
    colnames(fit$single_effect_lfsr) <- traits
    rownames(fit$alpha) <- paste0("L", 1:lmax)
    effects <- matrix(0, r, L)
    rownames(effects) <- traits
    colnames(effects) <- css

    table <- data.frame(matrix(as.numeric(NA),
        prod(length(conditions) * length(markers)), 8))
    colnames(table) <- c("trait", "variant", "pos",
        "effect", "z", "lfsr", "cs", "sentinel")
    table$trait <- rep(conditions, length(markers))
    table$variant <- rep(markers, each = length(conditions))
    table$pos <- rep(pos, each = length(conditions))
    table$sentinel <- 0
    for (i in 1:L){
        l <- css[i]
        j <- fit$sets$cs[[l]]
        b <- fit$b1_rescaled[l, j, ]
        if (conditional_effect){
            b <- b / fit$alpha[l, j]
        }
        variant_names <- markers[j]
        variant_idx <- which(table$variant %in% variant_names)
        table[variant_idx, "cs"] <- l
        table[variant_idx, "lfsr"] <-
          rep(fit$single_effect_lfsr[l,], length(variant_names))
        table[variant_idx, "effect"] <- c(t(b))
        if(!is.null(fit$z)){
            table[variant_idx, "z"] <- c(t(fit$z[j, ]))
        }
        max_idx <- which.max(fit$alpha[l, j])
        table[which(table$variant == variant_names[max_idx]), "sentinel"] <- 1
        effects[, i] <- ifelse(is.null(nrow(b)), b, b[max_idx, ])
    }
    table <- table[which(!is.na(table$cs)), ]
    if (!missing(poslim)) {
        rows1 <- which(table$pos >= poslim[1] & table$pos <= poslim[2])
        table <- table[rows1, ]
    }
    table$marker <- paste0(table$variant, "(", table$cs, ")")
    pdat_sentinel <- table[which(table$sentinel == 1), ]
    pdat_sentinel <- unique(pdat_sentinel[, c("variant", "pos", "marker")])
    pdat_sentinel$pip <- fit$pip[match(pdat_sentinel$variant,
        fit$variable_names)]
    if (sentinel_only){
        table <- table[which(table$sentinel == 1), ]
    }
    if (!missing(cs_plot))
      table <- table[which(table$cs %in% cs_plot),]
    table$cs    <- factor(table$cs)
    table$trait <- factor(table$trait,traits)
    
    # Remove from the effects plot any effects that don't meet the lfsr
    # cutoff.
    rows <- which(table$lfsr < lfsr_cutoff)
    table <- table[rows,]

    # Create the PIP plot.
    pip_plot <- ggplot(pdat, aes(x = "pos", y = "pip")) +
      geom_point(color = "darkblue", shape = 20, size = 1.25) +
      geom_point(shape = 1, size = 1.25, stroke = 1.25, data = pdat_cs,
                 mapping = aes(x = "pos", y = "pip", color = "cs")) +
      geom_text_repel(data = pdat_sentinel,
                      mapping = aes(x = "pos", y = "pip", label = "marker"),
                      size = 2.2, segment.size = 0.35, max.overlaps = Inf,
                      min.segment.length = 0) +
      xlim(poslim[1], poslim[2]) +
      scale_color_manual(values = cs_colors) +
      guides(color = guide_legend(override.aes = list(shape = 20,
                                                      size = 1.5))) +
      labs(x = sprintf("chromosome %d position (Mb)", chr),
           y = "PIP", color = "CS") +
      theme_cowplot(font_size = 9)

    # Create the effect plot.
    if (nrow(table) > 0) {
        table$effect_sign <- factor(table$effect > 0)
        table$effect_size <- abs(table$effect)
        levels(table$effect_sign) <- c("-1", "+1")
        effect_plot <- ggplot(table,
            aes(x = "marker", y = "trait", fill = "effect_sign",
                size = "effect_size")) +
            geom_point(shape = 21, stroke = 0.5, color = "white") +
            scale_y_discrete(drop = FALSE) +
            scale_fill_manual(values = c("darkblue", "red"), drop = FALSE) +
            scale_size(range = c(1, 5),
                breaks = unname(quantile(table$effect_size,
                    seq(0, 1, length.out = 4)))) +
            labs(x = "", y = "", fill = "effect sign", size = "effect size") +
            guides(fill = guide_legend(override.aes = list(size = 2)),
                size = guide_legend(override.aes = list(shape = 21, fill = "black",
                    color = "white", stroke = 0.5))) +
            theme_cowplot(font_size = 9) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                panel.grid = element_line(color = "lightgray", linetype = "dotted",
                    size = 0.3))
        # If requested, add colored dots to the top of the plot showing CS
        # membership.
        if (add_cs) {
            table$one <- 1
            p_cs <- ggplot(table, aes(x = "marker", y = "one", color = "cs")) +
                geom_point(shape = 20, size = 2.5) +
                scale_x_discrete(drop = FALSE) +
                scale_color_manual(values = cs_colors) +
                labs(x = "", y = "") +
                theme_cowplot(font_size = 12) +
                theme(axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    axis.line = element_blank())
            effect_plot <- plot_grid(p_cs, effect_plot, nrow = 2, ncol = 1,
                rel_heights = c(1, 2), axis = "lr", align = "v")
        }
    }else{
        effect_plot <- NULL
    }

    if (nrow(table) > 0 && all(!is.na(table$z))){
        table$z_sign <- factor(table$z > 0)
        table$z_size <- abs(table$z)
        levels(table$z_sign) <- c("-1", "+1")
        z_plot <- ggplot(table, aes(x = "marker", y = "trait",
            fill = "z_sign", size = "z_size")) +
            geom_point(shape = 21, stroke = 0.5, color = "white") +
            scale_y_discrete(drop = FALSE) +
            scale_fill_manual(values = c("darkblue", "red"), drop = FALSE) +
            scale_size(range = c(1, 5),
                breaks = unname(quantile(table$z_size,
                    seq(0, 1, length.out = 4)))) +
            labs(x = "", y = "", fill = "z-score sign", size = "z-score size") +
            guides(fill = guide_legend(override.aes = list(size = 2)),
                size = guide_legend(override.aes = list(shape = 21, fill = "black",
                    color = "white", stroke = 0.5))) +
            theme_cowplot(font_size = 9) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                panel.grid = element_line(color = "lightgray", linetype = "dotted",
                    size = 0.3))
        if (add_cs){
            z_plot <- plot_grid(p_cs, z_plot, nrow = 2, ncol = 1,
                rel_heights = c(1, 2), axis = "lr", align = "v")
        }
    }else{
        z_plot <- NULL
    }

    # Output the (1) PIP plot, (2) effect plot, (3) z scores plot and
    # (4) the table of effect estimates.
    return(list(pip_plot = pip_plot,
                effect_plot = effect_plot,
                z_plot = z_plot,
                effects = effects))
}
