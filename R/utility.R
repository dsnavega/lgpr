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
#' @noRd
#'
discrete_gaussian <- function(location, scale, interval, n) {

  # Prediction Domain
  domain <- seq(from = interval[1], to = interval[2], length.out = n)

  # Compute Discrete Normalization Factor (z)
  constant <- (1 / (scale * sqrt(2 * pi)))
  kernel <- sum(exp(-0.5 * ((domain - location) / scale) ^ 2))
  norm <- constant * kernel

  # Compute Discrete Probability Value
  constant <- (1 / (scale * sqrt(2 * pi) * norm))
  kernel <- exp(-0.5 * ((domain - location) / scale) ^ 2)
  probability <-  constant * kernel

  return(list(x = domain, y = probability))

}

#' @author David Senhora Navega
#' @noRd
#'
predict_discrete_gaussian <- function(object, alpha = 0.05) {

  x <- object$x
  y <- object$y

  z <- order(y / sum(y), decreasing = TRUE)
  area <- 0.0
  i <- 1

  while (i <= length(y)) {
    area <- sum(y[z[1:i]])
    if (area >= (1 - alpha))
      break

    i <- i + 1
  }

  x_max <- x[which.max(y)]
  x_range <- range(x[sort(z[1:i])])
  location <- c(x_max, x_range)
  names(location) <- c("estimate", "lower", "upper")

  return(location)

}

#' @author David Senhora Navega
#' @noRd
#'
local_gaussian <- function(location, scale, interval, n, alpha, digits = 3) {
  m <- length(location)
  local.list <- lapply(seq_len(m), function(i) {

    # Local Uncertainty using Discrete Gaussian
    local.gaussian <- discrete_gaussian(
      location = location[i], scale = scale[i],
      interval = interval, n = n
    )

    if (is.na(scale[i]))
      local.gaussian$y <- rep(0, times = length(local.gaussian$x))

    # Local Prediction
    if (is.na(location[i]) | is.na(scale[i])) {
      prediction <- rep(NA, 3)
    } else {
      prediction <- predict_discrete_gaussian(
        object = local.gaussian, alpha = alpha
      )
      prediction[1] <- location[i]
    }

    gaussian.list <- list(
      prediction = round(x = prediction, digits = digits),
      location = location[i], scale = scale[i],
      plot = local.gaussian, alpha = alpha
    )

    return(gaussian.list)

  })

  return(invisible(local.list))

}

#' @author David Senhora Navega
#' @noRd
#'
gaussian_prediction_matrix <- function(x) {
  return(do.call(rbind, lapply(x, function(x) x$prediction)))
}

#' Plot prediction object of memmento network when uncertainty is modelled from
#' a local gaussian around the point of interest.
#'
#' @author David Senhora Navega
#' @import ggplot2 stats
#'
#' @param object a list with the gaussian prediction of the memmento network
#' @param normalize a logical stating if the density should be normalized by its
#' maximum value. Default = T
#' @param digits an integer passed to the round() function. Default = 3.
#' @param label a character of the label of the x-axis.
#'
#'
#'

#' @author David Senhora Navega
#' @noRd
#'
plot_lgpr <- function(object, normalize = T, digits = 3, label = NULL) {

  if (is.na(object$location))
    return(NULL)

  # Helpers ----
  nearest_value <- function(x, z, ties = FALSE) {
    .nearest <- Vectorize(function(x) {
      if (ties) {
        which(abs(z - x) == min(abs(z - x)))
      } else {
        which.min(abs(x - z))
      }
    })
    return(.nearest(x))
  }

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

  lolli.index <- nearest_value(object$prediction, subset_df$x)
  lollipop_df <- data.frame(
    x = subset_df$x[lolli.index], y = subset_df$y[lolli.index]
  )

  .density <- stats::approxfun(
    x = x, y = if (normalize) { y / max(y) } else { y }
  )

  lollipop_df <- data.frame(
    x = x_values, y = .density(x_values)
  )

  grob <- ggplot2::ggplot(data = df, mapping = aes(x = x, y = y)) +
    ggplot2::geom_line(size = 0.75) +
    ggplot2::geom_area(data = subset_df, alpha = 0.25, linetype = "dashed") +
    ggplot2::geom_point(data = lollipop_df, size = 5.5, shape = 19) +
    ggplot2::geom_segment(
      data = lollipop_df, linetype = 2, size = 1.25,
      mapping = aes(x = x, xend = x, y = 0, yend = y)
    ) +
    ggplot2::scale_x_continuous(
      name = label,
      limits = range(x),
      breaks = seq(min(x) - (min(x) %% 10), max(x), 10),
      minor_breaks = seq(min(x) - (min(x) %% 5), max(x), 5),
    ) +
    ggplot2::ylab(label = y_label) +
    ggplot2::ggtitle(
      label = "Predictive Distribution (Local Gaussian)",
      subtitle = paste0(
        "Predicted: ",x_values[1]," [", x_values[2]," - ",x_values[3],"]\n",
        "Conditional Variance: ", round(object$scale, digits = digits),"\n",
        "Confidence: ", 1 - object$alpha
      )
    ) +
    ggplot2::theme_classic(base_line_size = 0.75) +
    ggplot2::theme(
      text = ggplot2::element_text(size = 14),
      axis.ticks = ggplot2::element_line(
        size = 1, colour = "black", lineend = "round",linetype = 3
      ),
      axis.title = ggplot2::element_text(family = "sans", face = "plain", size = 14),
      axis.text = ggplot2::element_text(
        family = "sans",face = "plain", size = 12, colour = "black"
      )
    )

  return(grob)

}
