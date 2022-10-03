# code to prepare objects to run occupancy models on high-performance computing

# load libraries
library(R2jags)
library(rslurm)
library(sparta)
library(reshape2)
library(coda)

# cleaned and tidied occurrence record data
# fake data example
taxa_data <- readRDS("data/example_taxa_data.rds")

visitData <- sparta::formatOccData(taxa = taxa_data$CONCEPT,
                                   site = taxa_data$TO_GRIDREF,
                                   survey = taxa_data$TO_STARTDATE)

# Set region aggregates
region_aggs <- readRDS("data/region_aggs_non-crossed.rds")
reg_data <- readRDS("data/regions_non-crossed.rds")

# Define function that loops through species
slurm_occDetFunc <- function(taxa_name){
  
  out <- sparta::occDetFunc(taxa_name = as.character(taxa_name),
                            occDetdata = visitData$occDetdata,
                            spp_vis = visitData$spp_vis,
                            write_results = TRUE,
                            n_chains = 3,
                            n_iterations = 32000,
                            burnin = 30000,
                            thinning = 6,
                            nyr = 2,
                            modeltype = c('ranwalk', 'halfcauchy', 'catlistlength'),
                            additional.parameters = "a",
                            regional_codes = reg_data,
                            region_aggs = region_aggs,
                            rem_aggs_with_missing_regions = FALSE,
                            allowSitesMultiRegions = TRUE,
                            provenance = "example data",
                            return_data = FALSE,
                            seed = 123)
  return(NULL)
}

# Create roster
pars <- data.frame(taxa_name = as.character(names(visitData[['spp_vis']])[-1]))

# Create the job script and the R script needed to run the process on 
# lotus using slurm. Note: you can edit the templates used. These are
# found in the slurm folder in your R library (run '.Library' to find).
# You will need to add the command to load jaspy: module add jaspy
sjob <- rslurm::slurm_apply(f = slurm_occDetFunc,
                            params = pars, 
                            jobname = 'example_SS',
                            nodes = nrow(pars), 
                            cpus_per_node = 1, 
                            submit = TRUE,
                            global_objects = c('visitData', 'reg_data', 'region_aggs'),
                            slurm_options = list(time = '23:59:00', 
                                                 mem = 8 * 1024,
                                                 partition = 'short-serial',
                                                 error = '%a.err'))
