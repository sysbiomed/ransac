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
                        family = 'binomial.glm', name = '',
                        baseline = list(), show.title = T,
                        only_consensus = T, outliers = NULL,
                        show.misclass = T,
                        print.plots = T,
                        ...) {
  strOut <- textConnection("foo", open = 'w', local = TRUE)
  sink(strOut)
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
  cat(sprintf('Results for k = %d iterations', length(result.ransac$error.array)), file = strOut, sep = '\n')

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
  if (FALSE && !only_consensus) {
    #
    if (any(!is.na(result.ransac$models$inliers.inliers.ydata))) {
      description$rans.5 <- 'RANSAC Initial + Inliers'
      models$rans.5      <- result.ransac$models$inliers.inliers.ydata
    }
    if (any(!is.na(result.ransac$models$inliers.all.ydata))) {
      description$rans.2 <- 'RANSAC Initial + {In,Out}liers'
      models$rans.2      <- result.ransac$models$inliers.all.ydata
    }
    if (any(!is.na(result.ransac$models$all.inliers.model.inliers))) {
      description$rans.3 <- 'RANSAC Refitted + Inliers'
      models$rans.3      <- result.ransac$models$all.inliers.model.inliers
    }
    if (any(!is.na(result.ransac$models$all.inliers.all.ydata))) {
      description$rans.4 <- 'RANSAC Refitted + {In,Out}liers'
      models$rans.4      <- result.ransac$models$all.inliers.all.ydata
    }
    if (any(!is.na(result.ransac$models$all.inliers.consensus))) {
      description$ransac <- 'RANSAC Refitted + Consensus'
    }
  } else {
    if (any(!is.na(result.ransac$models$all.inliers.consensus))) {
      description$ransac <- 'RANSAC'
    }
  }
  if (any(!is.na(result.ransac$models$all.inliers.consensus))) {
    models$ransac      <- result.ransac$models$all.inliers.consensus
  }

  # check if argument exists or is empty list
  if(is.null(baseline) || (is.list(baseline) && length(baseline) == 0)) {
    models$baseline      <- family.fun$fit.model(xdata, ydata, ...)
    description$baseline <- 'Baseline'
    if (any(models$baseline[[family]]$lambda == Inf)) {
      cat(sprintf('Error when fitting baseline'), file = strOut, sep = '\n')
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
    cat(sprintf('Using given baseline model'), file = strOut, sep = '\n')

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
  g1 <- ggplot(data = my.df.class %>% arrange(desc(type))) +
    geom_vline(aes(xintercept = (sum(ydata == 0) + 0.5)), color = 'black', linetype = 'dotted') +
    geom_text(aes(x=sum(ydata == 0)/2), y = 1.05, label= "Class 0", color = '#999999', family = 'Helvetica-Narrow', size = 3.5) +
    geom_text(aes(x=sum(ydata == 0) + sum(ydata == 1) / 2), y = 1.05, label= "Class 1", color = '#999999', family = 'Helvetica-Narrow', size = 3.5) +
    scale_y_continuous(expand = c(0,.1)) +
    theme_minimal() + #theme(legend.position="none") +
    ylab('Classification') +
    xlab('Model') +
    geom_point(aes(ix, value, color = type), alpha = 1, shape = 5) #+ facet_wrap( ~ type , ncol = 3)
  #
  if (show.title)      { g1 <- g1 + ggtitle(name) }
  if (FALSE && !only_consensus) { g1 <- g1 + theme(axis.ticks = element_blank(), axis.text.x = element_blank())}

  #
  # Error plot
  g2 <- ggplot(data = my.df) +
    theme_minimal() + theme(legend.position="none") +
    ylab('Error') + xlab('Model') +
    geom_quasirandom(aes(type, value, color = type), bandwidth = 2, method = 'quasirandom') +
    scale_y_continuous(limits = c(-1,1))
  if (show.title)      { g2 <- g2 + ggtitle(name) }
  if (FALSE && !only_consensus) { g2 <- g2 + theme(axis.ticks = element_blank(), axis.text.x = element_blank())}

  if (print.plots) {
    print(g1)
    print(g2)
  }

  #
  # Textual results
  #

  #
  cat(sprintf('Information on RANSAC and Baseline model'), file = strOut, sep = '\n')
  non.zero <- list()
  has.shown.sep <- F
  for(ix in results.names) {
    non.zero[[ix]] <- family.fun$coef(models[[ix]], ...)[-1] != 0
    if (ix %in% names(baseline)) {
      has.shown.sep <- T
      cat(sprintf('----------------- Baseline ----------'), file = strOut, sep = '\n')
    }
    cat(sprintf('  %d Co-Var. % 4d Obs %e %s in %s',
              sum(non.zero[[ix]]),
              family.fun$nobs(models[[ix]]),
              my.rmse[[ix]],
              family.fun$model.error.type,
              description[[ix]]), file = strOut, sep = '\n')
  }

  #
  # False Positive/Negative textual output
  #

  #
  cat(sprintf(''), file = strOut, sep = '\n')
  cat(sprintf(''), file = strOut, sep = '\n')
  cat(sprintf('False Positive/Negative'), file = strOut, sep = '\n')
  for(ix in results.names) {
    cat(sprintf('  % 3d / % 3d (total: % 3d) -- %s',
              sum(miscl[[ix]]$false.pos),
              sum(miscl[[ix]]$false.neg),
              sum(miscl[[ix]]$false.pos) + sum(miscl[[ix]]$false.neg),
              description[[ix]]), file = strOut, sep = '\n')
  }

  #
  outlier.other.name <- 'Other'
  if (is.null(outliers)) {
    outliers <- factor(array('', length(ydata)))
    outlier.other.name <- 'Miscl.'
  }
  cat(sprintf(''), file = strOut, sep = '\n')
  if (show.misclass) {
    cat(sprintf(''), file = strOut, sep = '\n')
    cat(sprintf('Misclassifications index'), file = strOut, sep = '\n')
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
      cat(sprintf(''), file = strOut, sep = '\n')
      cat(sprintf('  %s', description[[ix]]), file = strOut, sep = '\n')
      cat(sprintf('    False Positives % 4d', sum(miscl[[ix]]$false.pos)), file = strOut, sep = '\n')
      if (length(levels(outliers)) > 1) {
        cat(sprintf('           Natural (% 4d) %s', sum(natural.pos),   paste(which(natural.pos), collapse = ', ')), file = strOut, sep = '\n')
        cat(sprintf('           Pertur. (% 4d) %s', sum(perturbed.pos), paste(which(perturbed.pos), collapse = ', ')), file = strOut, sep = '\n')
        cat(sprintf('      Nat. + Pert. (% 4d) %s', sum(both.pos),      paste(which(both.pos), collapse = ', ')), file = strOut, sep = '\n')
      }
      cat(sprintf('             %s (% 4d) %s', outlier.other.name, sum(other.pos),     paste(which(other.pos), collapse = ', ')), file = strOut, sep = '\n')
      cat(sprintf('    False Negatives % 4d', sum(miscl[[ix]]$false.neg)), file = strOut, sep = '\n')
      if (length(levels(outliers)) > 1) {
        cat(sprintf('           Natural (% 4d) %s', sum(natural.neg),   paste(which(natural.neg), collapse = ', ')), file = strOut, sep = '\n')
        cat(sprintf('           Pertur. (% 4d) %s', sum(perturbed.neg), paste(which(perturbed.neg), collapse = ', ')), file = strOut, sep = '\n')
        cat(sprintf('      Nat. + Pert. (% 4d) %s', sum(both.neg),      paste(which(both.neg), collapse = ', ')), file = strOut, sep = '\n')
      }
      cat(sprintf('             %s (% 4d) %s', outlier.other.name, sum(other.neg),     paste(which(other.neg), collapse = ', ')), file = strOut, sep = '\n')
      }
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
  coef.df[is.na(coef.df)] <- 0
  cat(sprintf('Coefficients:'), file = strOut, sep = '\n')
  cat(capture.output(coef.df[,colSums(coef.df) != 0]), file = strOut, sep = '\n')

  sink()
  close(strOut)
  #
  # return the misclassifications
  return(list(debug = foo, description = description, misclassifications = miscl, plots = list(classification = g1, error = g2)))
}

