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
#'
#' @return A \code{ggplot} object.
#'
#' @export
#'
#' @examples
#' P <- matrix(c(0.5,0.3,0.8, 0.5,0.7,0.2), 3, 2)
#' plot_structure(P)
plot_structure <- function (P, pops = NULL,
                            label = NULL, map.indiv = NULL, map.pop = NULL, gap = NULL,
                            colors = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00",
                                       "#ffff33","#a65628","#f781bf","#999999"),
                            font.size = 9)
{
  if (is.null(pops))
  {
    pops <- order(colMeans(P))
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
  ggplot(dat,aes_string(x = "sample",y = "prop",color = "population",
                        fill = "population")) +
    geom_col() +
    scale_x_continuous(limits = c(0,max(dat$sample) + 1),breaks = ticks,
                       labels = names(ticks)) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    labs(x = "",y = "population proportion") +
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
                    population = rep(pops, each = n),
                    prop = c(P[, pops]))
  dat$population <- factor(dat$population, pops)
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
