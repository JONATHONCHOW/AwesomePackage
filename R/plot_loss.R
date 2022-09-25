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
#' @param L A list where each element is a vector of loss functions with respect to the number of iterations.
#' @param methods A vector of the same length as `L` where each element is a name of the fitting algorithm, such as "em", "sqp", "vi", "svi".
#' @param sample.rate The sampling rate of the loss function.
#' @param epsilon A small, positive number added to the vertical axis so that the
#'     logarithmic scale does not over-emphasize very small differences.
#' @param colors The colors used to draw loss curves.
#' @param linetypes The line types used to draw loss curves.
#' @param linesizes The line sizes used to draw loss curves.
#' @param shapes The shapes used to draw points at iterations.
#' @param fills The fill colors used to draw points at iterations.
#' @param theme The ggplot2 theme.
#' @param title Title of the plot, such as "K = 2", "K = 3".
#'
#' @return A \code{ggplot} object.
#'
#' @export
#'
#' @examples
#' L <- list(c(-100,-20,-10,-1,-0.5,-0.3,-0.1), c(-30,-5,-1,-0.1,-0.05))
#' plot_loss(L, c("fun","more fun"), 1)
plot_loss <- function (L, methods,
                       sample.rate, epsilon = 0.01,
                       colors = c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
                                  "#D55E00","#CC79A7"),
                       linetypes = "solid", linesizes = 0.5, shapes = 19, fills = "white",
                       theme = function() theme_cowplot(12),
                       title = NULL)
{
  n <- length(L)
  if (length(colors) < n)
    colors <- rep(colors,length.out = n)
  if (length(linetypes) < n)
    linetypes <- rep(linetypes,length.out = n)
  if (length(linesizes) < n)
    linesizes <- rep(linesizes,length.out = n)
  if (length(shapes) < n)
    shapes <- rep(shapes,length.out = n)
  if (length(fills) < n)
    fills <- rep(fills,length.out = n)
  ylab <- "distance from best loss function"
  xlab <- "iteration"
  dat <- compile_loss_data(L, methods, sample.rate, epsilon)
  ggplot(dat,aes_string(x = "iter",y = "dist.loss",color = "method",
                        linetype = "method",size = "method")) +
    geom_line(na.rm = TRUE) +
    geom_point(data = dat,
               mapping = aes_string(x = "iter",y = "dist.loss",color = "method",
                                    fill = "method",shape = "method"),
               inherit.aes = FALSE,na.rm = TRUE) +
    scale_y_continuous(trans = "log10") +
    scale_color_manual(values = colors) +
    scale_linetype_manual(values = linetypes) +
    scale_size_manual(values = linesizes) +
    scale_shape_manual(values = shapes) +
    scale_fill_manual(values = fills) +
    labs(x = xlab,y = ylab,title = title) +
    theme()
}

# Create a data frame suitable for loss plot.
compile_loss_data <- function (L, methods, sample.rate, epsilon)
{
  n <- length(L)
  dat <- list()
  for (i in 1:n)
  {
    dat[[i]] <- data.frame(iter = sample.rate * rep(1:length(L[[i]])),
                           dist.loss = L[[i]],
                           method = methods[i])
  }
  dat <- do.call(rbind,dat)
  dat$dist.loss <- max(dat$dist.loss) - dat$dist.loss + epsilon
  return(dat)
}

#' @title Plot the relationship between index and K
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
#' @description Draw a diagram of the relationship between the index and K using package ggplot2. The index can be loglikelihood, error, ELBO, and time.
#'
#' @param L A list where each element is a vector of index with respect to K.
#' @param methods A vector of the same length as `L` where each element is a name of the fitting algorithm, such as "em", "sqp", "vi", "svi".
#' @param index.id Choose index. Should be one of \code{"loglik"}, \code{"error"}, \code{"elbo"}, \code{"time"}.
#' @param start.point The initial point of K.
#' @param colors The colors used to draw curves.
#' @param linetypes The line types used to draw curves.
#' @param linesizes The line sizes used to draw curves.
#' @param shapes The shapes used to draw points.
#' @param fills The fill colors used to draw points.
#' @param theme The ggplot2 theme.
#' @param title Title of the plot.
#'
#' @return A \code{ggplot} object.
#'
#' @export
#'
#' @examples
#' L <- list(c(-50,-20,-10,-1,-0.5), c(-30,-5,-1,-0.1))
#' plot_index_vs_K(L, c("fun","more fun"), index.id = "loglik")
#' L <- list(c(0.1,0.5,0.7,1.2), c(1.6,0.2,0.8,1.5,2.4))
#' plot_index_vs_K(L, c("fun","more fun"), index.id = "error")
#' L <- list(c(-10,-2,-0.5,-0.1), c(-5,-1,-0.1))
#' plot_index_vs_K(L, c("fun","more fun"), index.id = "elbo")
#' L <- list(c(10,15,20), c(12,18,30,32))
#' plot_index_vs_K(L, c("fun","more fun"), index.id = "time")
plot_index_vs_K <- function (L, methods, index.id = c("loglik","error","elbo","time"),
                             start.point = 2,
                             colors = c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
                                        "#D55E00","#CC79A7"),
                             linetypes = "solid", linesizes = 0.5, shapes = 19, fills = "white",
                             theme = function() theme_cowplot(12),
                             title = NULL)
{
  n <- length(L)
  if (length(colors) < n)
    colors <- rep(colors,length.out = n)
  if (length(linetypes) < n)
    linetypes <- rep(linetypes,length.out = n)
  if (length(linesizes) < n)
    linesizes <- rep(linesizes,length.out = n)
  if (length(shapes) < n)
    shapes <- rep(shapes,length.out = n)
  if (length(fills) < n)
    fills <- rep(fills,length.out = n)
  data <- compile_index_vs_K_data(L, methods, index.id, start.point)
  dat <- data$dat
  ylab <- data$ylab
  xlab <- "K"
  ggplot(dat,aes_string(x = "K",y = "index",color = "method",
                        linetype = "method",size = "method")) +
    geom_line(na.rm = TRUE) +
    geom_point(data = dat,
               mapping = aes_string(x = "K",y = "index",color = "method",
                                    fill = "method",shape = "method"),
               inherit.aes = FALSE,na.rm = TRUE) +
    scale_color_manual(values = colors) +
    scale_linetype_manual(values = linetypes) +
    scale_size_manual(values = linesizes) +
    scale_shape_manual(values = shapes) +
    scale_fill_manual(values = fills) +
    labs(x = xlab,y = ylab,title = title) +
    theme()
}

# Create a data frame suitable for index vs K plot.
compile_index_vs_K_data <- function (L, methods, index.id, start.point)
{
  n <- length(L)
  dat <- list()
  if (index.id == "loglik")
  {
    ylab <- "loglikelihood"
  }
  if (index.id == "error")
  {
    ylab <- "error"
  }
  if (index.id == "elbo")
  {
    ylab <- "ELBO"
  }
  if (index.id == "time")
  {
    ylab <- "time to reach convergence"
  }
  for (i in 1:n)
  {
    dat[[i]] <- data.frame(K = rep(start.point:(length(L[[i]]) + start.point - 1)),
                           index = L[[i]],
                           method = methods[i])
  }
  dat <- do.call(rbind,dat)
  return(list(dat = dat, ylab = ylab))
}
