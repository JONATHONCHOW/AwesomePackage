#' @title Plot the population proportion
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_col
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_text
#' @importFrom cowplot theme_cowplot
#'
#' @description Plot the population proportion of individuals using package ggplot2.
#'
#' @param P The proportion matrix.
#' @param pops Population order options.
#' @param label The original order of individuals. This option is only for data that needs to be grouped.
#' @param map.indiv The new order of individuals. This option is only for data that needs to be grouped.
#' @param map.pop The order of populations. This option is only for data that needs to be grouped.
#' @param gap Gaps between groups. This option is only for data that needs to be grouped.
#' @param colors Theme color options.
#' @param font.size Font size used in plot.
#' @param title Title of the plot, such as "EM", "SQP", "VI", "SVI".
#' @param subtitle Subtitle of the plot, such as "K = 2", "K = 3".
#'
#' @return A \code{ggplot} object.
#'
#' @export
#'
#' @examples
#' P <- matrix(c(0.5,0.3,0.8, 0.5,0.7,0.2), 3, 2)
#' plot_structure(P, title = "FUN")
plot_structure <- function (P, pops = NULL,
                            label = NULL, map.indiv = NULL, map.pop = NULL, gap = NULL,
                            colors = c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
                                       "#EE2C2C","#CC79A7","#8968CD","#FF83FA","#EECFA1",
                                       "#A52A2A","#4169E1","#FFFF00","#BFEFFF","#FF1493"),
                            font.size = 9,
                            title = NULL, subtitle = NULL)
{
  if (is.null(pops))
  {
    pops <- order(colMeans(P), decreasing = T)
  }
  if (is.null(label))
  {
    dat <- compile_structure_data(P, pops)
    ticks <- NULL
  }
  else
  {
    group <- compile_group_structure_data(P, pops, label, map.indiv, map.pop, gap)
    dat <- group$dat
    ticks <- group$ticks
  }
  ggplot(dat,aes_string(x = "sample",y = "prop",color = "pop",
                        fill = "pop")) +
    geom_col() +
    scale_x_continuous(limits = c(0,max(dat$sample) + 1),breaks = ticks,
                       labels = names(ticks)) +
    scale_y_continuous(breaks = NULL) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    labs(x = NULL,y = NULL,title = title,subtitle = subtitle) +
    theme_cowplot(font.size) +
    theme(axis.line   = element_blank(),
          axis.ticks  = element_blank(),
          axis.text.x = element_text(angle = 45,hjust = 1))
}

# Create a data frame suitable for structure plot.
compile_structure_data <- function (P, pops)
{
  n <- nrow(P)
  k <- length(pops)
  dat <- data.frame(sample = rep(1:n, times = k),
                    pop = rep(pops, each = n),
                    prop = c(P[, pops]))
  dat$pop <- factor(dat$pop, pops)
  return(dat)
}

# Create a grouped data frame suitable for structure plot.
compile_group_structure_data <- function (P, pops, label, map.indiv, map.pop, gap)
{
  Q <- matrix(nrow = nrow(P), ncol = ncol(P))
  for (i in 1:length(map.indiv))
  {
    Q[i, ] <- P[which(label == map.indiv[i]), ]
  }
  groups <- unique(map.pop)
  ticks <- rep(0, length(groups))
  names(ticks) <- groups
  dat <- NULL
  m <- 0
  for (j in groups) {
    i          <- which(map.pop == j)
    out        <- compile_structure_data(Q[i,], pops)
    out$sample <- out$sample + m
    n          <- length(i)
    dat        <- rbind(dat, out)
    ticks[j]   <- m + n/2
    m          <- m + n + gap
  }
  return(list(dat = dat,ticks = ticks))
}
