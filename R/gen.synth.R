gen.synth <- function(obs,
                      new.coef,
                      perturbation.perct = 0.1,
                      perfect.model = F,
                      xdata.fun = rnorm,
                      perturbation.random = FALSE) {

  len <- obs

  # intercept
  intercept <- new.coef[1]
  new.coef  <- new.coef[-1]

  # get number of variables from coefficient
  nvars     <- length(new.coef)

  # generate a random X with normal distribution
  xdata.all <- sapply(seq(nvars), function(e, len) { xdata.fun(len)}, len)

  # http://stats.stackexchange.com/questions/124538/how-to-generate-a-large-full-rank-random-correlation-matrix-with-some-strong-cor#125020
  #k <- 1
  #sigma.m <- matrix(rnorm(nvars * k), nrow = nvars, ncol = k)
  #S = sigma.m %*% t(sigma.m) + diag(runif(nvars))
  #S = diag(1/sqrt(diag(S))) %*% S %*% diag(1/sqrt(diag(S)));

  #out <- xdata.all %*% chol(S)
  #cor(xdata.all)
  #hist(cor(out)[lower.tri(cor(out), 0)], breaks = 50)

  # sigma.m <- matrix(runif(nvars^2, min = 0, max = 1) * 10,nrow=nvars)
  # diag(sigma.m) <- array(1, nvars)
  # random.ix <- sample(which(upper.tri(sigma.m, 0)), ceiling(nvars * .1))
  # sigma.m[random.ix]
  # sigma.m[random.ix] <- sigma.m[random.ix] * 10
  # sigma.m[random.ix]
  # sigma.m[lower.tri(sigma.m, 0)] <- t(sigma.m)[lower.tri(sigma.m, 0)]

  #xdata.all <- exp(rmvnorm(100,mean=rep(0, nvars), sigma= sigma.m, method="svd"))

  # add names
  rownames(xdata.all) <- seq(nrow(xdata.all))
  colnames(xdata.all) <- paste0('X', seq(ncol(xdata.all)))

  # linear combination with a bias
  z  <- intercept + xdata.all %*% new.coef

  # pass through an inv-logit function
  pr <- 1 / (1 + exp(-z))

  # bernoulli response variable
  ydata  <- data.frame(prob = pr, real_class = rbinom(len, 1 , pr))
  xdata <- xdata.all
  if (perfect.model) {
    ydata$real_class = (ydata$prob > 0.5) * 1
  }

  # rename rows
  rownames(ydata) <- rownames(xdata) <- seq(obs)

  # perturbation (balanced)

  # start with real class
  ydata$logit_class <- ydata$real_class
  # calculate size of perturbation
  pert.size <- round(obs * perturbation.perct)
  if (perturbation.random) {
    # TRUE: class 1
    # FALSE: class 0
    start.by <- (runif(1) > 0.5)
    my.sample <- list(class1 = sample(which(ydata$real_class == 1), sum((ydata$real_class == 1))),
                      class0 = sample(which(ydata$real_class == 0), sum((ydata$real_class == 0))))
    for (ix in seq(pert.size)) {
      ix.name <- sprintf('class%d', start.by * 1)
      ydata$logit_class[my.sample[[ix.name]][1]] <- (ydata$logit_class[my.sample[[ix.name]][1]] - 1) * -1
      #
      my.sample[[ix.name]] <- my.sample[[ix.name]][-1]
      start.by <- !start.by
    }
  } else {
    pr.ix <- sort(pr,index.return=T)$ix

    tail.size <- ceiling(obs * perturbation.perct / 2)

    if (pert.size > 0) {
      forced.outliers.count <- c(0,0)
      for (ix in seq(obs)) {
        first.ix <- pr.ix[ix]
        last.ix  <- pr.ix[obs - (ix - 1)]
        first.pr <- pr[first.ix]
        last.pr  <- pr[last.ix]
        # skip if natural outlier
        #  or already made perturbations on class 0
        if (first.pr - ydata$real_class[first.ix] < 0.5 &&
            forced.outliers.count[1] <= tail.size) {
          ydata$logit_class[first.ix] <- (ydata$real_class[first.ix] - 1) * -1
          forced.outliers.count[1] <- forced.outliers.count[1] + 1
        }
        # skip if natural outlier
        #  or already made perturbations on class 1
        if (last.pr - ydata$real_class[last.ix] < 0.5 &&
            forced.outliers.count[2] <= tail.size) {
          ydata$logit_class[last.ix] <- (ydata$real_class[last.ix] - 1) * -1
          forced.outliers.count[2] <- forced.outliers.count[2] + 1

        }
        #
        if (sum(forced.outliers.count) >= pert.size) {
          break
        }
      }
    }
  }

  ydata$outliers = array('inlier', length(ydata$real_class))
  # natural outliers
  ydata$outliers[ydata$real_class != (ydata$prob > 0.5) * 1] <- 'natural'
  ydata$outliers[ydata$real_class != ydata$logit_class] <- 'perturbation'

  # return a list of xdata and ydata
  return(list(xdata = xdata, ydata = ydata))
}
