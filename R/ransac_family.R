#' Get RANSAC family function
#'
#' @param family
#'
#' @return a ransac.family function
#' @export
#'
#' @examples
#' ransac.family('ransac.binomial.glm')
ransac.family <- function(family) {
  if (is.character(family)) {
    family.fun <- switch(family,
                         binomial.glmnet        = ransac.binomial.glmnet(auc = F, residuals = 'pearson'),
                         binomial.glmnet.se     = ransac.binomial.glmnet(auc = F, residuals = 'squared.error'),
                         binomial.glmnet.auc    = ransac.binomial.glmnet(auc = T, residuals = 'pearson'),
                         binomial.glmnet.se.auc = ransac.binomial.glmnet(auc = T, residuals = 'squared.error'),
                         #
                         binomial.glm        = ransac.binomial.glm(auc = F, residuals = 'pearson'),
                         binomial.glm.se     = ransac.binomial.glm(auc = F, residuals = 'squared.error'),
                         binomial.glm.auc    = ransac.binomial.glm(auc = T, residuals = 'pearson'),
                         binomial.glm.se.auc = ransac.binomial.glm(auc = T, residuals = 'squared.error'))
  } else {
    family.fun <- family
  }
  if (is.null(family.fun)) {
    stop('family is not well defined, see documentation')
  }
  return(family.fun)
}
