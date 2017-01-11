print.ransac.results <- function(results.ransac, xdata, ydata, family = 'binomial', name = '') {
  #
  # retrieve family by name
  if (is.character(family)) {
    family.fun <- switch(family,
                         binomial = family.binomial())
  } else {
    family.fun <- family
  }
  if (is.null(family.fun)) {
    stop('family is not well defined, see documentation')
  }
  #
  if (results.ransac$best.error == Inf) {
    flog.info('Could not find a good fit')
  } else {
    models <- list()
    models$ransac <- results.ransac$best.model
    models$glmnet <- family.fun$fit.model(xdata, ydata)
    #
    my.fits <- list()
    my.fits$ransac <- family.fun$predict(models$ransac, xdata)
    my.fits$glmnet <- family.fun$predict(models$glmnet, xdata)
    #
    my.rmse <- list()
    my.rmse$ransac <- family.fun$model.error(my.fits$ransac, ydata)
    my.rmse$glmnet <- family.fun$model.error(my.fits$glmnet, ydata)
    #
    my.df <- data.frame()
    my.df <- rbind(my.df, data.frame(ix = seq(nrow(xdata)), val = (ydata - my.fits$ransac), type = rep('ransac', nrow(xdata))))
    my.df <- rbind(my.df, data.frame(ix = seq(nrow(xdata)), val = (ydata - my.fits$glmnet), type = rep('glmnet', nrow(xdata))))
    my.df$type <- factor(my.df$type)
    colnames(my.df) <- c('ix', 'value', 'type')
    #
    my.df.class <- data.frame()
    ix.sorted <- sort(ydata, index.return = T)$ix
    my.df.class <- rbind(my.df.class, data.frame(ix = seq(nrow(xdata)), val = (my.fits$ransac[ix.sorted]), type = rep('ransac', nrow(xdata))))
    my.df.class <- rbind(my.df.class, data.frame(ix = seq(nrow(xdata)), val = (my.fits$glmnet[ix.sorted]), type = rep('glmnet', nrow(xdata))))
    my.df.class <- rbind(my.df.class, data.frame(ix = seq(nrow(xdata)), val = (ydata[ix.sorted]),          type = rep('ydata',  nrow(xdata))))
    my.df.class$type <- factor(my.df.class$type)
    colnames(my.df.class) <- c('ix', 'value', 'type')
    #
    # misclassifications
    miscl <- list()
    miscl$ransac <- list()
    miscl$glmnet <- list()
    miscl$ransac$false.pos <- sum(round(ydata - my.fits$ransac) == -1)
    miscl$ransac$false.neg <- sum(round(ydata - my.fits$ransac) == 1)
    #
    miscl$glmnet$false.pos <- sum(round(ydata - my.fits$glmnet) == -1)
    miscl$glmnet$false.neg <- sum(round(ydata - my.fits$glmnet) == 1)
    #
    title.str <- sprintf('%s : ransac FP %d / %d FN | glmnet FP %d / %d FN',
                         name,
                         miscl$ransac$false.pos,
                         miscl$ransac$false.neg,
                         miscl$glmnet$false.pos,
                         miscl$glmnet$false.neg)
    #
    g <- ggplot(data = my.df.class) + theme_minimal() + ggtitle(title.str) + ylab('Classification') + xlab('Model')
    g <- g + geom_point(aes(ix, value, color = type))
    print(g)
    #
    g <- ggplot(data = my.df) + theme_minimal() + ggtitle(title.str) + ylab('Error') + xlab('Model')
    g <- g + geom_quasirandom(aes(type, value, color = type))
    print(g)
    #
    flog.info('RMSE:\n   ransac = %f\n   glmnet = %f', my.rmse$ransac, my.rmse$glmnet)
  }
}
