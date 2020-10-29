# This file is part of lgpr
#
# <lgpr: regression uncertainty from local gaussian distributions>
# Copyright (C) 2020, David Senhora Navega
#
# lgpr is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# lgpr is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with lgpr. If not, see <http://www.gnu.org/licenses/>.
#
# David Senhora Navega
# Laboratory of Forensic Anthropology
# Department of Life Sciences
# University of Coimbra
# Calçada Martim de Freitas, 3000-456, Coimbra
# Portugal

#' @author David Senhora Navega
#' @import truncnorm
#' @noRd
#'
truncate_gaussian <- function(location, scale, interval, n, alpha) {

  if (is.na(location) | is.na(scale)) {
    prediction <- rep(NA, 3)
    x.axis <- seq(from = interval[1], to = interval[2], length.out = n)
    y.axis <- rep(x = 0, times = n)
    object <- list(
      prediction = prediction,
      location = location, scale = scale,
      plot = list(x = x.axis, y = y.axis),
      alpha = alpha
    )
    return(object)
  }

  # Confidence
  p.alpha <- c(alpha / 2, 1 - (alpha / 2))

  # Prediction Interval
  bound <- truncnorm::qtruncnorm(
    p = p.alpha, a = interval[1], b = interval[2], mean = location, sd = scale
  )

  # Prediction Domain
  x.axis <- seq(from = interval[1], to = interval[2], length.out = n)
  x.axis <- sort(c(bound, x.axis))

  # Gaussian Density
  y.axis <- truncnorm::dtruncnorm(
    x = x.axis, a = interval[1], b = interval[2], mean = location, sd = scale
  )

  # Prediction
  prediction <- c(x.axis[which.max(y.axis)], bound)
  names(prediction) <- c("estimate", "lower", "upper")

  # Object
  object <- list(
    prediction = prediction,
    location = location, scale = scale,
    plot = list(x = x.axis, y = y.axis),
    alpha = alpha
  )

  return(object)

}

#' @author David Senhora Navega
#' @noRd
#'
local_gaussian <- function(location, scale, interval, n, alpha, digits = 3) {

  local.list <- lapply(seq_len(length.out = length(location)), function(ith) {
    gaussian.list <- truncate_gaussian(
      location = location[ith], scale = scale[ith],
      interval = interval, n = n, alpha = alpha
    )
    return(gaussian.list)
  })

  return(invisible(local.list))

}

#' @author David Senhora Navega
#' @noRd
#'
lgpr_prediction_matrix <- function(x) {
  return(do.call(rbind, lapply(x, function(x) x$prediction)))
}

#' @author David Senhora Navega
#' @noRd
#'
plot_lgpr <- function(object, normalize = T, digits = 3, label = NULL) {

  if (is.na(object$location))
    return(NULL)

  x <- round(object$plot$x, digits = digits)
  y <- object$plot$y

  if (normalize) {
    df <- data.frame(x = c(min(x),x,max(x)), y = c(0, y / max(y), 0))
    y_label <- "Density (Normalized)"
  } else {
    df <- data.frame(x = c(min(x),x,max(x)), y = c(0, y, 0))
    y_label <- "Density"
  }

  x_values <- round(object$prediction, digits = digits)

  subset_df <- df[df$x >= x_values[2] & df$x <= x_values[3],]

  lollipop_df <- data.frame(
    x = subset_df$x[subset_df$x %in% x_values],
    y = subset_df$y[subset_df$x %in% x_values]
  )

  x_values <- format(x_values, nsmall = digits)
  p <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_line(size = 0.75) +
    ggplot2::geom_area(data = subset_df, alpha = 0.25, linetype = "dashed") +
    ggplot2::geom_point(data = lollipop_df, size = 5, shape = 19) +
    ggplot2::geom_segment(
      data = lollipop_df, linetype = 2, size = 1,
      mapping = ggplot2::aes(x = x, xend = x, y = 0, yend = y)
    ) +
    ggplot2::scale_x_continuous(
      name = label,
      limits = range(x),
      breaks = seq(min(x) - (min(x) %% 10), max(x), 10),
      minor_breaks = seq(min(x) - (min(x) %% 5), max(x), 5),
    ) +
    ggplot2::ylab(label = y_label) +
    ggplot2::ggtitle(
      label = "Predictive Distribution",
      subtitle = paste0(
        "Predicted: ",x_values[1]," [", x_values[2]," - ",x_values[3],"]\n",
        "Conditional Variance: ", round(object$scale, digits = digits),"\n",
        "Confidence: ", 1 - object$alpha
      )
    ) +
    ggplot2::theme_classic(base_line_size = 0.75) +
    ggplot2::theme(
      text = ggplot2::element_text(size = 13),
      axis.ticks = ggplot2::element_line(
        size = 1, colour = "black", lineend = "round",linetype = 3
      ),
      axis.title = ggplot2::element_text(family = "sans", face = "plain", size = 14),
      axis.text = ggplot2::element_text(
        family = "sans",face = "plain", size = 11, colour = "black"
      )
    )

  return(p)

}
