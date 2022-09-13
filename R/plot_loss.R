#' @title Plot the loss function
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_linetype_manual
#' @importFrom ggplot2 scale_size_manual
#' @importFrom ggplot2 scale_shape_manual
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 labs
#' @importFrom cowplot theme_cowplot
#'
#' @description Plot the loss function against the number of iterations using package ggplot2.
#'
#' @param L A vector represents the loss function.
#' @param sample.rate The sampling rate of the loss function.
#' @param method The name of the fitting algorithm, such as "em", "sqp", "vi", "svi".
#' @param epsilon A small, positive number added to the vertical axis so that the
#'     logarithmic scale does not over-emphasize very small differences.
#' @param colors The colors used to draw loss curves.
#' @param linetypes The line types used to draw loss curves.
#' @param linesizes The line sizes used to draw loss curves.
#' @param shapes The shapes used to draw points at the selected iterations.
#' @param fills The fill colors used to draw points at the selected iterations.
#' @param theme The \sQuote{ggplot2} \dQuote{theme}.
#'
#' @return A \code{ggplot} object.
#'
#' @export
#'
#' @examples
#' L <- c(-100,-20,-10,-1,-0.5,-0.3,-0.1)
#' plot_loss(L, 1, "fun")
plot_loss <- function (L, sample.rate, method,
                       epsilon = 0.01,
                       colors = c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
                                  "#D55E00","#CC79A7"),
                       linetypes = "solid", linesizes = 0.5, shapes = 19, fills = "white",
                       theme = function() theme_cowplot(12))
{
  dat <- compile_loss_data(L, sample.rate, epsilon)
  xlab <- "iteration"
  ylab <- "distance from best loss function"
  ggplot(dat,aes_string(x = dat$iter,y = dat$dist.loss,color = "method",
                                linetype = "method",size = "method")) +
           geom_line(na.rm = TRUE) +
           geom_point(data = dat,
                      mapping = aes_string(x = dat$iter,y = dat$dist,color = "method",
                                           fill = "method",shape = "method"),
                      inherit.aes = FALSE,na.rm = TRUE) +
           scale_y_continuous(trans = "log10") +
           scale_color_manual(values = colors) +
           scale_linetype_manual(values = linetypes) +
           scale_size_manual(values = linesizes) +
           scale_shape_manual(values = shapes) +
           scale_fill_manual(values = fills) +
           labs(x = xlab,y = ylab) +
           theme()
}

# Create a data frame suitable for loss plot.
compile_loss_data <- function (L, sample.rate, epsilon)
{
  dat <- data.frame(iter = sample.rate * rep(0:(length(L) - 1)),
                    dist.loss = max(L) - L + epsilon)
  return(dat)
}
