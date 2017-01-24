#' Function to print the results from RANSAC
#'
#' @param result.ransac
#' @param xdata
#' @param ydata
#' @param family
#' @param name
#' @param ydata.original
#' @param show.title
#' @param only_consensus
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot.ransac <- function(result.ransac, xdata, ydata,
                        ydata.prob = NULL,
                        family = 'binomial', name = '',
                        baseline = list(), show.title = T,
                        only_consensus = T, outliers = NULL,
                        ...) {
  #
  # retrieve family by name
  if (is.character(family)) {
    family.fun <- ransac.family(family)
  } else {
    family.fun <- family
  }
  if (is.null(family.fun)) {
    stop('family is not well defined, see documentation')
  }

  # probabilities vector
  #  defaults to 0/1 if not given
  if (is.null(ydata.prob)) {
    ydata.prob <- ydata
  }

  #
  flog.info('Results for k = %d iterations', length(result.ransac$error.array))

  # Populate penalty.factor as 1 if is not defined
  if (!exists('penalty.factor')) {
    penalty.factor <- array(1, ncol(xdata))
  }

  #
  # Build models
  #  - if no baseline model is passed in arguments the
  #     a simple fit with xdata / ydata is performed
  #  - otherwise, build models from baseline and result.ransac

  description <- list()
  models <- list()

  # build from result.ransac
  if (!only_consensus) {
    #
    description$rans.5 <- 'RANSAC Initial + Inliers'
    models$rans.5      <- result.ransac$models$inliers.inliers.ydata
    description$rans.2 <- 'RANSAC Initial + {In,Out}liers'
    models$rans.2      <- result.ransac$models$inliers.all.ydata
    description$rans.3 <- 'RANSAC Refitted + Inliers'
    models$rans.3      <- result.ransac$models$all.inliers.model.inliers
    description$rans.4 <- 'RANSAC Refitted + {In,Out}liers'
    models$rans.4      <- result.ransac$models$all.inliers.all.ydata
    description$ransac <- 'RANSAC Refitted + Consensus'
  } else {
    description$ransac <- 'RANSAC'
  }
  models$ransac      <- result.ransac$models$all.inliers.consensus

  # check if argument exists or is empty list
  if(is.null(baseline) || (is.list(baseline) && length(baseline) == 0)) {
    models$baseline      <- family.fun$fit.model(xdata, ydata, ...)
    description$baseline <- 'Baseline'
    if (any(models$baseline[[family]]$lambda == Inf)) {
      flog.info('Error when fitting baseline')
      return()
    }
  }
  # a model is given
  else if (is.list(baseline) && !('list' %in% class(baseline))) {
    models$baseline      <- baseline
    description$baseline <- 'Baseline'
  }
  # else build from baseline argument
  else {
    flog.info('Using given baseline model')

    for(ix in names(baseline)) {
      models[[ix]]      <- baseline[[ix]]
      #
      description[[ix]] <- gsub('[._-]', ' ', ix)
    }
  }

  #
  # build list of names that will show in the plots and results
  #  this is done by joining the baseline model(s) with the output
  #  from ransac
  results.names <- names(models)

  # fit the xdata to all the models
  my.fits <- list()
  for(ix in results.names) {
    my.fits[[ix]] <- family.fun$predict(models[[ix]], xdata, ...)
  }

  # find model error
  #  TODO: change name from my.rmse to my.model.error
  my.rmse <- list()
  for(ix in results.names) {
    my.rmse[[ix]] <- family.fun$model.error(my.fits[[ix]], ydata)
  }
  names(my.rmse) <- results.names

  # build the error plot's data
  my.df <- data.frame()
  for(ix in results.names) {
    my.df <- rbind(my.df,
                   data.frame(ix = seq(nrow(xdata)),
                              val = (ydata.prob - my.fits[[ix]]),
                              type = rep(description[[ix]], nrow(xdata))))
  }
  my.df$type <- factor(my.df$type)
  colnames(my.df) <- c('ix', 'value', 'type')

  # build the classification plot's data
  my.df.class <- data.frame()
  ix.sorted <- sort(ydata, index.return = T)$ix
  for(ix in results.names) {
    my.df.class <- rbind(my.df.class,
                         data.frame(ix = seq(nrow(xdata)),
                                    val = (my.fits[[ix]][ix.sorted]),
                                    type = rep(description[[ix]], nrow(xdata))))
  }
  # add the observed values
  my.df.class <- rbind(my.df.class,
                       data.frame(ix = seq(nrow(xdata)),
                                  val = (ydata[ix.sorted]),
                                  type = rep('Observed',  nrow(xdata))))
  my.df.class$type <- factor(my.df.class$type)
  colnames(my.df.class) <- c('ix', 'value', 'type')

  #
  # calculate misclassifications
  miscl <- list()
  for(ix in results.names) {
    miscl[[ix]] <- list()
    miscl[[ix]]$false.pos <- round(ydata - my.fits[[ix]]) == -1
    miscl[[ix]]$false.neg <- round(ydata - my.fits[[ix]]) == 1
  }

  #
  # Plot the results
  #

  #
  # Classification plot
  g <- ggplot(data = my.df.class) +
    theme_minimal() +
    ylab('Classification') +
    xlab('Model') +
    geom_point(aes(ix, value, color = type)) + facet_wrap( ~ type , ncol = 2)
  #
  if (show.title)      { g <- g + ggtitle(name) }
  if (!only_consensus) { g <- g + theme(axis.ticks = element_blank(), axis.text.x = element_blank())}
  print(g)

  #
  # Error plot
  g <- ggplot(data = my.df) +
    theme_minimal() +
    ylab('Error') + xlab('Model') +
    geom_quasirandom(aes(type, value, color = type), bandwidth = 2, method = 'quasirandom') +
    scale_y_continuous(limits = c(-1,1))
  if (show.title)      { g <- g + ggtitle(name) }
  if (!only_consensus) { g <- g + theme(axis.ticks = element_blank(), axis.text.x = element_blank())}
  print(g)

  #
  # Textual results
  #

  #
  flog.info('Information on RANSAC and Baseline model')
  non.zero <- list()
  has.shown.sep <- F
  for(ix in results.names) {
    non.zero[[ix]] <- family.fun$coef(models[[ix]], ...)[-1] != 0
    if (ix %in% names(baseline)) {
      has.shown.sep <- T
      flog.info('----------------- Baseline ----------')
    }
    flog.info('  %d Co-Var. % 4d Obs %e %s in %s',
              sum(non.zero[[ix]]),
              family.fun$nobs(models[[ix]]),
              my.rmse[[ix]],
              family.fun$model.error.type,
              description[[ix]])
  }

  #
  # False Positive/Negative textual output
  #

  #
  flog.info('')
  flog.info('')
  flog.info('False Positive/Negative')
  for(ix in results.names) {
    flog.info('  % 3d / % 3d (total: % 3d) -- %s',
              sum(miscl[[ix]]$false.pos),
              sum(miscl[[ix]]$false.neg),
              sum(miscl[[ix]]$false.pos) + sum(miscl[[ix]]$false.neg),
              description[[ix]])
  }

  #
  if (is.null(outliers)) {
    outliers <- factor(array('', length(ydata)))
  }
  flog.info('')
  flog.info('')
  flog.info('Misclassifications index')
  for(ix in results.names) {
    natural.pos <- miscl[[ix]]$false.pos & outliers == 'Natural'
    natural.neg <- miscl[[ix]]$false.neg & outliers == 'Natural'
    #
    perturbed.pos <- miscl[[ix]]$false.pos & outliers == 'Perturbed'
    perturbed.neg <- miscl[[ix]]$false.neg & outliers == 'Perturbed'
    #
    both.pos <- miscl[[ix]]$false.pos & outliers == 'Natural and Perturbed'
    both.neg <- miscl[[ix]]$false.neg & outliers == 'Natural and Perturbed'
    #
    other.pos <- miscl[[ix]]$false.pos & outliers != 'Natural and Perturbed' & outliers != 'Natural' & outliers != 'Perturbed'
    other.neg <- miscl[[ix]]$false.neg & outliers != 'Natural and Perturbed' & outliers != 'Natural' & outliers != 'Perturbed'
    #
    {
    flog.info('')
    flog.info('  %s', description[[ix]])
    flog.info('    False Positives % 4d', sum(miscl[[ix]]$false.pos))
    flog.info('           Natural (% 4d) %s', sum(natural.pos),   paste(which(natural.pos), collapse = ', '))
    flog.info('           Pertur. (% 4d) %s', sum(perturbed.pos), paste(which(perturbed.pos), collapse = ', '))
    flog.info('      Nat. + Pert. (% 4d) %s', sum(both.pos),      paste(which(both.pos), collapse = ', '))
    flog.info('             Other (% 4d) %s', sum(other.pos),     paste(which(other.pos), collapse = ', '))
    flog.info('    False Negatives % 4d', sum(miscl[[ix]]$false.pos))
    flog.info('           Natural (% 4d) %s', sum(natural.neg),   paste(which(natural.neg), collapse = ', '))
    flog.info('           Pertur. (% 4d) %s', sum(perturbed.neg), paste(which(perturbed.neg), collapse = ', '))
    flog.info('      Nat. + Pert. (% 4d) %s', sum(both.neg),      paste(which(both.neg), collapse = ', '))
    flog.info('             Other (% 4d) %s', sum(other.neg),     paste(which(other.neg), collapse = ', '))
    }
  }

  #
  # build a data.frame with coefficients
  coef.df <- data.frame()
  for (ix in names(models)) {
    coef.df <- rbind(coef.df, as.vector(family.fun$coef(models[[ix]], ...)))
  }
  rownames(coef.df) <- unlist(description)
  colnames(coef.df) <- c('Intercept', colnames(xdata))
  flog.info('Coefficients:', coef.df, capture = T)

  #
  # return the misclassifications
  return(list(description = description, misclassifications = miscl))
}

