dataset.brca.tnbc <- function(select = c(), glmnet.alpha = NULL, mc.cores = 1, convert.to.1.0 = T) {
  #
  data("fpkm.per.tissue", 'clinical', package = 'brca.data')
  #

  # select columns of clinical data to find triple negative biomarker
  tnbc.vars <- c('breast_carcinoma_estrogen_receptor_status',
                 'breast_carcinoma_progesterone_receptor_status',
                 'lab_proc_her2_neu_immunohistochemistry_receptor_status')
  clinical.tnbc <- clinical$primary.solid.tumor[, tnbc.vars]

  # replace all '' or NA with 'Indeterminate'
  clinical.tnbc[is.na(clinical.tnbc) | clinical.tnbc == '' | clinical.tnbc == 'Equivocal'] <- 'Indeterminate'

  #
  # Find all that have at least one positve and keep all
  #  After find all with Indeterminate and exclude from dataset

  # With at least one positve are marked as not TNBC
  tnbc.status.pos <- clinical.tnbc == 'Positive'
  not.tnbc.id     <- rownames(tnbc.status.pos[which(sapply(seq(nrow(tnbc.status.pos)), function(ix) { any(tnbc.status.pos[ix,]) })),])

  # With at least one Indeterminate
  tnbc.status.ind <- clinical.tnbc[!(rownames(clinical.tnbc) %in% not.tnbc.id),] == 'Indeterminate'
  not.tnbc.id     <- rownames(tnbc.status.ind[which(sapply(seq(nrow(tnbc.status.ind)), function(ix) { any(tnbc.status.ind[ix,]) })),])
  clinical.tnbc <- clinical.tnbc[!(rownames(clinical.tnbc) %in% not.tnbc.id), ]

  tnbc.status.neg <- clinical.tnbc == 'Negative'
  tnbc.ix <- sapply(seq(nrow(clinical.tnbc)), function(ix) { all(tnbc.status.neg[ix,]) })

  # set cases with TNBC
  clinical.tnbc$tnbc <- 'NO_TNBC'
  clinical.tnbc[tnbc.ix, 'tnbc'] <- 'TNBC'

  # build ydata
  ydata <- data.frame(tnbc = clinical.tnbc$tnbc, row.names = rownames(clinical.tnbc))
  ydata$tnbc <- factor(ydata$tnbc)
  #
  xdata <- t(fpkm.per.tissue$primary.solid.tumor[,-which(duplicated(fpkm.per.tissue.barcode$primary.solid.tumor))])
  common.rows <- getParticipantCode(rownames(xdata)) %in% rownames(ydata)
  xdata <- xdata[common.rows,]
  #
  xdata.sd <- sapply(seq(ncol(xdata)), function(ix) {sd(xdata[,ix])})
  xdata <- xdata[,xdata.sd != 0]

  # normalize by
  for(ix in seq(ncol(xdata))) {
    xdata[,ix] <- (xdata[,ix] - mean(xdata[,ix])) / sd(xdata[,ix])
  }
  #
  if (length(select) > 0) {
    xdata <- xdata[, select]
    flog.info('Selecting only %d co-variates from \'select\' argument.')
  } else if (!is.null(glmnet.alpha)) {
    # balanced folds
    class.tnbc.ix <- which(ydata == 'TNBC')
    class.norm.ix <- which(ydata != 'TNBC')
    #
    foldid.struct <- balanced.cv.folds(seq_along(class.tnbc.ix), seq_along(class.norm.ix), nfolds = 10)
    foldid <- array(0, nrow(ydata))
    foldid[class.tnbc.ix] <- foldid.struct$output[[1]]
    foldid[class.norm.ix] <- foldid.struct$output[[2]]
    #
    my.fit <- cv.glmnet(xdata, as.vector(ydata == 'TNBC') * 1, family = 'binomial',
                        alpha = glmnet.alpha, nlambda = 100,
                        foldid = foldid,
                        mc.cores = mc.cores)
    flog.info('cv.glmnet selected %d variables', sum(coef(my.fit, s = 'lambda.min') != 0))
    selected.ix <- as.vector(coef(my.fit, s = 'lambda.min') != 0)
    # gotta take intercept
    selected.ix <- which(selected.ix)[-1] - 1
    xdata <- xdata[,selected.ix]
  } else {
    # do nothing
  }
  #
  if (convert.to.1.0) {
    ydata <- (ydata == 'TNBC') * 1
  }
  colnames(ydata) <- 'logit_class'
  #
  return(list(xdata = xdata, ydata = ydata, name = 'Breast Invasive Carcinome TNBC (triple negative breast cancer)'))
}
