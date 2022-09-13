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
#' @param colors Theme color options.
#' @param ticks Grouping options.
#' @param font.size Font size used in plot.
#'
#' @return A \code{ggplot} object.
#'
#' @export
#'
#' @examples
#' P <- matrix(c(0.5,0.3,0.8, 0.5,0.7,0.2), 3, 2)
#' pops <- order(colMeans(P))
#' colors <- c("red", "yellow")
#' plot_structure(P, pops, colors)
plot_structure <- function (P, pops, colors, ticks = NULL,
                            font.size = 9)
{
  dat <- compile_structure_data(P, pops)
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
  dat <- data.frame(sample = rep(1:n,times = k),
                    population = rep(pops,each = n),
                    prop = c(P[,pops]))
  dat$population <- factor(dat$population,pops)
  return(dat)
}
