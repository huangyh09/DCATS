
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
