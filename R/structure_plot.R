#' @title Plot the ancestral proportions
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
#' @description Plot the ancestral proportions of individuals using package ggplot2.
#'
#' @param L The proportion matrix.
#' @param topics Topic order options.
#' @param colors Theme color options.
#' @param ticks Grouping options.
#' @param font.size Font size used in plot.
#'
#' @export
#'
#' @examples
#' P <- matrix(c(0.5,0.5, 0.3,0.7, 0.8,0.2), 3, 2)
#' topics <- order(colMeans(P))
#' colors <- c("red", "yellow")
#' structure_plot(P, topics, colors)
structure_plot <- function (L, topics, colors, ticks = NULL,
                            font.size = 9)
{
  dat <- compile_structure_plot_data(L, topics)
  ggplot(dat,aes_string(x = "sample",y = "prop",color = "topic",
                        fill = "topic")) +
  geom_col() +
  scale_x_continuous(limits = c(0,max(dat$sample) + 1),breaks = ticks,
                     labels = names(ticks)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = "",y = "topic proportion") +
  theme_cowplot(font.size) +
  theme(axis.line   = element_blank(),
        axis.ticks  = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1))
}

# Create a data frame suitable for structure plot.
compile_structure_plot_data <- function (L, topics)
{
  n <- nrow(L)
  k <- length(topics)
  dat <- data.frame(sample = rep(1:n,times = k),
                    topic  = rep(topics,each = n),
                    prop   = c(L[,topics]))
  dat$topic <- factor(dat$topic,topics)
  return(dat)
}
