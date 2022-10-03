library(base, quietly = TRUE)
library(methods, quietly = TRUE)
library(datasets, quietly = TRUE)
library(utils, quietly = TRUE)
library(grDevices, quietly = TRUE)
library(graphics, quietly = TRUE)
library(stats, quietly = TRUE)
library(remotes, quietly = TRUE)
library(sf, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(coda, quietly = TRUE)
library(rjags, quietly = TRUE)
library(R2jags, quietly = TRUE)
library(sp, quietly = TRUE)
library(raster, quietly = TRUE)
library(HDInterval, quietly = TRUE)
library(effsize, quietly = TRUE)
library(betapart, quietly = TRUE)
library(adespatial, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(cowplot, quietly = TRUE)
library(wrappeR, quietly = TRUE)
library(gghalves, quietly = TRUE)
library(occAssess, quietly = TRUE)
library(rslurm, quietly = TRUE)
library(Matrix, quietly = TRUE)
library(lme4, quietly = TRUE)
library(sparta, quietly = TRUE)
library(reshape2, quietly = TRUE)

load('add_objects.RData')

.rslurm_func <- readRDS('f.RDS')
.rslurm_params <- readRDS('params.RDS')
.rslurm_more_args <- readRDS('more_args.RDS')
.rslurm_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
.rslurm_istart <- .rslurm_id * 1 + 1
.rslurm_iend <- min((.rslurm_id + 1) * 1, nrow(.rslurm_params))
.rslurm_result <- do.call(parallel::mcmapply, c(
    FUN = .rslurm_func,
    .rslurm_params[.rslurm_istart:.rslurm_iend, , drop = FALSE],
    MoreArgs = list(.rslurm_more_args),
    mc.cores = 1,
    mc.preschedule = TRUE,
    SIMPLIFY = FALSE))

saveRDS(.rslurm_result, file = paste0('results_', .rslurm_id, '.RDS'))
