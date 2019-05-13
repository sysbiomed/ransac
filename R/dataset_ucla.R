
#' Logistic dataset from UCLA
#'
#' from [UCLA](http://statistics.ats.ucla.edu/stat/r/dae/logit.htm)
#'
#' @return
#' @export
#'
#' @examples
dataset.ucla <- function() {
  # http://statistics.ats.ucla.edu/stat/r/dae/logit.htm
  mydata <- read.csv("http://www.ats.ucla.edu/stat/data/binary.csv")
  #
  xdata <- as.matrix(mydata[,c('gpa','rank','gre')])
  ydata <- data.frame(mydata$admit)
  colnames(ydata) <- 'logit_class'
  #
  return(list(xdata = xdata, ydata = ydata, name = 'UCLA'))
}
