#' Create a sequence of lambda values
#'
#' This generates a lambda vector with values that
#'  will be used in GLMNET in the coordinate descent
#'  algorightm. It is needed as it builds on each lambda
#'
#' @param lambda target that should be calculated
#' @param extend
#'
#' @return
#' @export
#'
#' @examples
#' find.lambda(0.08)
find.lambda <- function(lambda, extend = FALSE) {
  if (extend) {
    extend.vec <- seq(1, 3, 0.05) * lambda
  } else {
    extend.vect <- c()
  }
  ret.vect <- unique(sort(c(10, 8, 5,
                            c(seq(5, 50, 5), 4, 3, 2.5, 2, 1.75, 1.5, 1.25, 1,
                              extend.vect
                            ) * lambda), decreasing = T))
  return(ret.vect)
}
