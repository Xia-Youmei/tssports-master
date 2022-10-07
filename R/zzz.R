# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Add 100 to the given number
#'
#' @param x a number, float or int
#'
#' @return a number
#' @export
#'
#' @examples
#' add100(5)
#' add100(-90)
add100 <- function(x) {
  x + 100
}

