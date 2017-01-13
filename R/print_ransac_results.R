print.ransac.results <- function(results.ransac, xdata, ydata, family = 'binomial', name = '', show.title = T) {
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
    models$rans.2 <- results.ransac$best.model.all
    models$glmnet <- family.fun$fit.model(xdata, ydata, mc.cores = 14, nlambda = 100, alpha = 0.7)
    #
    my.fits <- list()
    my.fits$ransac <- family.fun$predict(models$ransac, xdata)
    my.fits$rans.2 <- family.fun$predict(models$rans.2, xdata)
    my.fits$glmnet <- family.fun$predict(models$glmnet, xdata)
    #
    my.rmse <- list()
    my.rmse$ransac <- family.fun$model.error(my.fits$ransac, ydata)
    my.rmse$rans.2 <- family.fun$model.error(my.fits$rans.2, ydata)
    my.rmse$glmnet <- family.fun$model.error(my.fits$glmnet, ydata)
    #
    my.df <- data.frame()
    my.df <- rbind(my.df, data.frame(ix = seq(nrow(xdata)), val = (ydata - my.fits$ransac), type = rep('RANSAC', nrow(xdata))))
    my.df <- rbind(my.df, data.frame(ix = seq(nrow(xdata)), val = (ydata - my.fits$rans.2), type = rep('RANSAC (All)', nrow(xdata))))
    my.df <- rbind(my.df, data.frame(ix = seq(nrow(xdata)), val = (ydata - my.fits$glmnet), type = rep('Baseline', nrow(xdata))))
    my.df$type <- factor(my.df$type)
    colnames(my.df) <- c('ix', 'value', 'type')
    #
    my.df.class <- data.frame()
    ix.sorted <- sort(ydata, index.return = T)$ix
    my.df.class <- rbind(my.df.class, data.frame(ix = seq(nrow(xdata)), val = (my.fits$ransac[ix.sorted]), type = rep('RANSAC', nrow(xdata))))
    my.df.class <- rbind(my.df.class, data.frame(ix = seq(nrow(xdata)), val = (my.fits$rans.2[ix.sorted]), type = rep('RANSAC (All)', nrow(xdata))))
    my.df.class <- rbind(my.df.class, data.frame(ix = seq(nrow(xdata)), val = (my.fits$glmnet[ix.sorted]), type = rep('Baseline', nrow(xdata))))
    my.df.class <- rbind(my.df.class, data.frame(ix = seq(nrow(xdata)), val = (ydata[ix.sorted]),          type = rep('Observed',  nrow(xdata))))
    my.df.class$type <- factor(my.df.class$type)
    colnames(my.df.class) <- c('ix', 'value', 'type')
    #
    # misclassifications
    miscl <- list()
    miscl$ransac <- list()
    miscl$rans.2 <- list()
    miscl$glmnet <- list()
    miscl$ransac$false.pos <- sum(round(ydata - my.fits$ransac) == -1)
    miscl$ransac$false.neg <- sum(round(ydata - my.fits$ransac) == 1)
    #
    miscl$rans.2$false.pos <- sum(round(ydata - my.fits$rans.2) == -1)
    miscl$rans.2$false.neg <- sum(round(ydata - my.fits$rans.2) == 1)
    #
    miscl$glmnet$false.pos <- sum(round(ydata - my.fits$glmnet) == -1)
    miscl$glmnet$false.neg <- sum(round(ydata - my.fits$glmnet) == 1)
    #
    title.str <- sprintf('%s : RANSAC FP %d / %d FN | RANSAC All FP %d / %d |\n  FN Ridge FP %d / %d FN',
                         name,
                         miscl$ransac$false.pos,
                         miscl$ransac$false.neg,
                         miscl$rans.2$false.pos,
                         miscl$rans.2$false.neg,
                         miscl$glmnet$false.pos,
                         miscl$glmnet$false.neg)
    #
    g <- ggplot(data = my.df.class) + theme_minimal() + ylab('Classification') + xlab('Model')
    if (show.title) { g <- g + ggtitle(title.str) }
    g <- g + geom_point(aes(ix, value, shape = type, color = type))
    print(g)
    #
    g <- ggplot(data = my.df) + theme_minimal() + ylab('Error') + xlab('Model')
    if (show.title) { g <- g + ggtitle(title.str) }
    g <- g + geom_quasirandom(aes(type, value, shape = type, color = type))
    print(g)
    #
    non.zero.ransac <- coef(models$ransac, s = 'lambda.min') != 0
    non.zero.rans.2 <- coef(models$rans.2, s = 'lambda.min') != 0
    non.zero.glmnet <- coef(models$glmnet, s = 'lambda.min') != 0
    count.common <- sum(rownames(non.zero.glmnet)[as.vector(non.zero.glmnet)] %in% rownames(non.zero.ransac)[as.vector(non.zero.ransac)])
    count.comm.2 <- sum(rownames(non.zero.glmnet)[as.vector(non.zero.glmnet)] %in% rownames(non.zero.rans.2)[as.vector(non.zero.rans.2)])
    count.comm.3 <- sum(rownames(non.zero.ransac)[as.vector(non.zero.ransac)] %in% rownames(non.zero.rans.2)[as.vector(non.zero.rans.2)])
    #
    flog.info(paste0('Information on RANSAC and Baseline model:\n',
                     '            Co-Var. in RANSAC: %d\n',
                     '      Co-Var. in RANSAC (All): %d\n',
                     '          Co-Var. in Baseline: %d\n',
                     '  Common co-var.     (G vs R): %d\n',
                     '  Common co-var. (G vs R.all): %d\n',
                     '  Common co-var. (R vs R.all): %d\n',
                     '----------------- RMSE --------------\n',
                     '                  RMSE RANSAC: %g\n',
                     '                RMSE Baseline: %g'),
              sum(non.zero.ransac),
              sum(non.zero.glmnet),
              sum(non.zero.rans.2),
              count.common,
              count.comm.2,
              count.comm.3,
              my.rmse$ransac, my.rmse$glmnet)
  }
}
