
#' Convert a matrix to data frame
#' @param X A matrix of values
#' @export
#' @example
#' mtx_to_df(matrix(seq(12), nrow=3))
mtx_to_df <- function(X) {
    if (is.null(row.names(X))) {
        row.names(X) <- seq(nrow(X))
    }
    if (is.null(colnames(X))) {
        colnames(X) <- seq(ncol(X))
    }
    data.frame(Var1 = rep(row.names(X), ncol(X)),
               Var2 = rep(colnames(X), each = nrow(X)),
               value = c(as.matrix(X)))
}


#' Plot heatmap from a matrix
#'
#' @param mat A matrix to show, column by x-axis and row by y-axis
#' @param base_size Numeric value for the base size in theme_bw
#' @param digits Integer value for the number of digits to show
#' @param show_value Logical value for showing the value for each element or not
#' @import ggplot2
#' @export
heat_matrix <- function(mat, base_size=12, digits=2, show_value=TRUE){
    df <- mtx_to_df(mat)
    if (!is.null(rownames(mat))) {
        df$Var1 <- factor(df$Var1, rownames(mat))
    }
    if (!is.null(colnames(mat))) {
        df$Var2 <- factor(df$Var2, colnames(mat))
    }

    df$value <- round(df$value, digits = digits)
    heat.plot <- ggplot(df, aes_string(x = "Var1", y = "Var2")) +
        geom_tile(aes_string(fill = "value"), colour = "grey") +
        scale_fill_gradient(low = "white", high = "steelblue") +
        #theme_grey(base_size = base_size) +
        theme_bw(base_size = base_size) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        theme(legend.position = "none",
              panel.grid.major = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              axis.ticks.x = ggplot2::element_blank(),
              axis.ticks.y = ggplot2::element_blank())
    if (show_value) {
        heat.plot <- heat.plot +
            geom_text(aes_string(label = "value"),
                      vjust = 0.5, size=base_size*0.25)
    }
    heat.plot
}


#' Vocalo plot differential abundance test on effect size and p values
#' @param res_df, data.frame as output from `dcats_fit()` or `betabinLRT()` with
#' variables `pvals`, `rownames`, `prop2_mean`, `coeff_mean`, and `coeff_std`
#' @param min_pval, A float to cape the minimum of p values
#' @param min_prop, A float to cape the minimum of cell type proportion in
#' control
#' @param show_EffectStd, A bool for indicating if showing the standard error of
#' effect size
#'
#' @import ggrepel, ggplot2
#' @export
#'
volcano_plot <- function(df_plot, min_pval=10^(-10), min_prop=3*10^(-4),
                         show_EffectStd=TRUE) {
    df_plot["cellType"] <- factor(rownames(df_plot))
    df_plot$pvals[df_plot$pvals < min_pval] = min_pval
    df_plot$prop2_mean[df_plot$prop2_mean < min_prop] = min_prop

    pp <- ggplot(df_plot, aes(x = coeff_mean, y = -log10(pvals),
                              color=cellType)) +
        geom_point(aes(color = cellType, size= prop2_mean)) +
        geom_hline(yintercept=-log10(0.05), color="firebrick") +
        geom_vline(xintercept=0) +
        ggrepel::geom_text_repel(aes(label = row.names(df_plot))) +
        guides(color=FALSE) + theme_bw() +
        labs(x="Effect size on logit", size="Prop_ctrl")
    if (show_EffectStd) {
        pp <- pp + geom_pointrange(aes(xmin=coeff_mean - coeff_std,
                                       xmax=coeff_mean + coeff_std))
    }

    pp
}


#' A dotplot for matrices with variable size and color
#'
#' @param size_mat A matrix of dot sizes with shape (n_row, n_column)
#' @param color_mat A matrix of dot colors with shape (n_row, n_column)
#' @param size_title A string value for the legend title on size
#' @param color_title A string value for the legend title on color
#' @param size_range A range value for size range in display
#' @param color_limits A vector or two floats for color limits in display. If
#' NULL, c(-1, 1) * max(abs(df_dat$Color)) is used.
#'
#' @return a ggplot object
#'
#' @export
#' @import ggplot2
#'
matrix_dotplot <- function(size_mat, color_mat,
                           size_title='Size', color_title='Color',
                           size_range=range(0.1, 4), color_limits=NULL) {
    # Check row and column names
    if (is.null(rownames(size_mat)))
        rownames(size_mat) <- paste0('Row', seq_len(nrow(size_mat)))
    if (is.null(colnames(size_mat)))
        colnames(size_mat) <- paste0('Column', seq_len(ncol(size_mat)))

    if (is.null(rownames(color_mat)))
        rownames(color_mat) <- paste0('Row', seq_len(nrow(color_mat)))
    if (is.null(colnames(color_mat)))
        colnames(color_mat) <- paste0('Column', seq_len(ncol(color_mat)))

    size_df = reshape2::melt(size_mat)
    color_df = reshape2::melt(color_mat)

    colnames(size_df) <- c('Row', 'Column', 'Size')
    colnames(color_df) <- c('Row', 'Column', 'Color')

    df_dat <- cbind(size_df, color_df['Color'])

    if (is.null(color_limits)) {
        color_limits = c(-1, 1) * max(abs(df_dat$Color))
    }

    df_dat$Color[df_dat$Color < color_limits[1]] <- color_limits[1]
    df_dat$Color[df_dat$Color > color_limits[2]] <- color_limits[2]

    ggplot(data = df_dat, aes(x = Column, y = Row)) +
        geom_point(aes(size = Size, colour = Color)) +
        xlab('') + ylab('') +
        theme_classic(base_line_size=0) +
        scale_size(range = size_range, name = size_title) +
        scale_color_distiller(palette = "RdBu", limits = color_limits,
                              name = color_title,
                              breaks=c(color_limits, color_limits/2, 0),
                              labels=c(paste('<', color_limits[1]),
                                       paste('>', color_limits[2]),
                                       color_limits/2, 0))
    # scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",
    #                       high = "red", space = "Lab", name = color_title)
}


#' Quantile-quantile plot for p value calibration
#'
#' @param pvalues P value vector or matrix
#' @return a ggplot object
#' @export
qq_plot <- function(pvalues) {
    obs_pval <- as.numeric(pvalues)
    exp_pval <- (rank(obs_pval, ties.method="first")+.5) / (length(obs_pval)+1)

    df.fig <- data.frame(obs_pval = obs_pval, exp_pval = exp_pval)

    ggplot(data = df.fig, aes(x = -log10(exp_pval), y = -log10(obs_pval))) +
        geom_point(shape = 1) +
        geom_abline(intercept = 0, slope = 1, color = "grey") +
        theme_bw()
}


