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

#' Local Gaussian Prediction for Regression Uncertainty Modelling
#'
#' @author David Senhora Navega
#' @export
#'
#' @param location a numeric vector for defining the location of the gaussian
#' distribution. See Details.
#' @param scale a numeric vector for defining the location of the gaussian
#' distribution. See Details.
#' @param interval a vector with two numeric values defining the upper and lower
#' bound of the discrete gaussian distribution.
#' @param n number of uniform steps defining the values at which the density of
#' the distribution (gaussian) is evaluated. Default is 512.
#' @param alpha tolerance of the confidence interval created by the gaussian
#' distribution. confidence = 1 - alpha.
#'
#' @returns a lgpr object
#'
#' @details
#' The location argument of the gaussian distribution in this context refers to
#' the point estimate or prediction obtained from a regression algorithm and the
#' scale argument refers to an estimate of the variance (error) associated to
#' such point estimate. If the regression algorithm do not output a variance
#' estimate by default, the variance can be model from the residuals of the model
#' for instance is up to the user / analyst to come up with a strategy to model
#' the variance (scale) associated with the regression predictions (location).
#' This function is intended to ease up the computation of prediction interval
#' when assuming a gaussian distribution for the residuals of the regression model
#' is reasonable. If the scale argument is composed of different values the
#' predictive intervals generated will be heteroscedastic. That is the reason
#' why this approach is named Local Gaussian Prediction because the predictive
#' interval are adjusted to difficulty of the prediction itself through the scale
#' argument.
#' In this implementation it is assumed that the prediction domain is bounded and
#' the limits are given by the interval argument. The prediction interval is
#' computed from the truncated gaussian distribution by a discretization
#' approximation.
#'
lgpr <- function(location, scale, interval, n = 512, alpha = 0.05)   {

  lgpr.object <- local_gaussian(
    location = location, scale = scale,
    interval = interval, alpha = alpha, n = n
  )

  object <- structure(.Data = lgpr.object, class = "lgpr")
  return(object)

}

#' Predict method for lgpr object
#'
#' @author David Senhora Navega
#' @export
#'
#' @param object an lgpr object
#' @param ... ....
#'
#' @return a 3-column matrix with prediction interval defined by the gaussian
#' uncertainty model.
#'
predict.lgpr <- function(object, ...) {

  if (!is.lgpr(object))
    stop("\n(-) x was not created by lpgr().")

  predictions <- gaussian_prediction_matrix(x = object)
  return(predictions)

}

#' Plot method for lgpr object
#' @author David Senhora Navega
#'
#' @export
#' @import ggplot2
#'
#' @param x a lgpr object
#' @param index the index of the prediction to be plotted. Default is 1.
#' @param normalize a logical stating if the density (y-axis) should be
#' normalized (values are divided by the maximum value). Default is TRUE
#' @param digits an integer defining the precision of the round function.
#' Default is 3.
#' @param label a character defining the label of the x-axis. Default is NULL.
#' @param ... ...
#'
#' @return a ggplot2 graphical object.
#'
plot.lgpr <- function(
  x, index = 1, normalize = T, digits = 3, label = NULL, ...
) {

  if (length(x) < index)
    stop("\n(-) lgpr model only contains ", length(x), " observation(s).")

  grob <- plot_lgpr(
    object = x[[index]], normalize = normalize, digits = digits, label = label
  )

  return(grob)

}

#' @author David Senhora Navega
#' @noRd
print.lgpr <- function(x) {
  cat("\n Regression Uncertainty Modelling using Gaussian Distribution")
  cat("\n Observations:", length(x))
}

#' @author David Senhora Navega
#' @noRd
is.lgpr <- function(object) {
  inherits(object, what = "lgpr")
}
