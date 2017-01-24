
library(futile.logger)

run.tnbc <- function(my.params,
                     rmd = file.path('vignettes', 'tnbc.Rmd'),
                     prefix = 'Report',
                     out.dir = 'output') {
  #'
  #' The goal of this script is to run RANSAC and save files
  devtools::load_all('.')

  #
  out.name <- sprintf('%s-%s-%d-%d-%d-%f-%f-%s',
                      prefix, my.params$data.origin, my.params$k, my.params$n, my.params$my.seed,
                      my.params$threshold, my.params$good.fit.perct, my.params$my.family)

  #' Create output directories
  if (!dir.exists(out.dir)) {
    dir.create(file.path(out.dir))
  }
  if (!dir.exists(file.path(out.dir, 'old'))) {
    dir.create(file.path(out.dir, 'old'))
  }

  #
  out.dir.this <- file.path(out.dir, prefix, out.name)
  if (dir.exists(out.dir.this)) {
    file.rename(out.dir.this,
                file.path(out.dir,
                          'old',
                          sprintf('%s-%05d-%s',
                                  format(Sys.time(), "%Y-%m-%d-%H:%M:%S"),
                                  round(runif(1) * 10000),
                                  out.name)))
    unlink(out.dir.this, recursive = T)
  }
  dir.create(out.dir.this, recursive = T)
  out.file <- file.path(out.dir.this, paste0(out.name, '.md'))

  flog.info('Directory: %s', out.dir.this)
  flog.info('File: %s', out.file)

  #'
  #' Call knit
  rmarkdown::render(rmd, 'all', output_file = out.file, output_dir = out.dir.this, params = my.params)
}
