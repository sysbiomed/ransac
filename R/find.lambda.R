#
#
#
find.lambda <- function(lambda) {
  return(unique(sort(c(1000,100, 10,
                c(50, 10, 5, 4, 3, 2, 1.5, 1,
                  #.9, .8, .7, .6, .5, .2, .3, .4, .1,
                  .5, .1,
                  1e-2, #5e-2,
                  1e-3, #5e-3,
                  1e-4, #5e-4,
                  1e-5, #5e-5,
                  1e-6, #5e-6,
                  1e-7#, 5e-7
                  ) * lambda), decreasing = T)))
}
