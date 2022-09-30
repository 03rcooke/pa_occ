## title: pa_occ
## full title: "Protected areas support more species than unprotected areas in Great Britain, but lose them equally rapidly"
## author: "Rob Cooke"
## date: "30/09/2022"

## Below is all the code used to prepare the data from raw and run the analyses, however the raw data is not available to share (owned by the data producer)
## So if you are interested in the analyses run, see the uncommented code below

#### Set-up ####
  
## Here we load the necessary packages

# # install JAGS in terminal
# sudo apt install jags

library(remotes) # remotes: package installation
library(sf) # sf: spatial manipulation
library(dplyr) # dplyr: data manipulation
library(tidyr) # tidyr: data manipulation
library(rjags) # rjags: occupancy models
library(R2jags) # R2jags: occupancy models
library(raster) # raster: spatial manipulation
library(HDInterval) # HDInterval: credible intervals
library(effsize) # effsize: effect sizes
library(betapart) # betapart: temporal beta diversity
library(adespatial) # adespatial: temporal beta diversity
library(ggplot2) # ggplot2: plotting
library(cowplot) # cowplot: plotting

# # sparta package from github
# withr::with_envvar(c(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true"),
# remotes::install_github("BiologicalRecordsCentre/sparta"))
# library(sparta) # sparta: occupancy models

# # wrappeR package from github
# remotes::install_github("https://github.com/BiologicalRecordsCentre/wrappeR", ref = "main")
library(wrappeR) # wrappeR: multi-species indicators

# # install occAssess from github
# remotes::install_github("https://github.com/robboyd/occAssess")
library(occAssess)

# # commit packages to renv library
# renv::snapshot()

# source functions
source("functions/pa_pr_fun.R") # source pa_pr function
source("functions/beta_tbi_func.R") # source beta_tbi_func function
source("functions/est_plot_func.R") # source est_plot_func function

# load_rdata function
# loads an RData file, and assigns it to an object name
load_rdata <- function(fileName) {
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# set location of temporary directory for raster package
rasterOptions(tmpdir = "temp")

# set plotting theme
ggplot2::theme_set(cowplot::theme_cowplot())

# global parameters
startyear <- 1990
endyear <- 2018
ci <- 0.95 # credible interval

#### Preprocessing ####

## Here we preprocess the protected areas data to create a raster of protection for GB

# # grid references
# # load: gr_ref
# gr_ref <- read.csv("data/gr_ref.csv") %>%
#   dplyr::mutate(easting = easting + 500,
#                 northing = northing + 500)

# # british national grid crs
# gbcrs <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs"
# 
# # protected areas (downloaded 13/10/2020 - https://www.protectedplanet.net/country/GBR)
# 
# # load protected area shapefiles
# pa0 <- sf::st_read("data/wdpa_0/WDPA_WDOECM_GBR_shp-polygons.shp") # too big to upload
# pa1 <- sf::st_read("data/wdpa_1/WDPA_WDOECM_GBR_shp-polygons.shp") # too big to upload
# pa2 <- sf::st_read("data/wdpa_2/WDPA_WDOECM_GBR_shp-polygons.shp") # too big to upload
# 
# pa <- rbind(pa0, pa1, pa2)
# 
# # load: pa
# pa <- pa %>%
#   # terrestrial only
#   dplyr::filter(MARINE == 0) %>%
#   # tier 1 sites only
#   dplyr::filter(DESIG_ENG %in% c("National Nature Reserve",
#                              "Ramsar Site, Wetland of International Importance",
#                              "Local Nature Reserve",
#                              "Site Of Special Scientific Interest (Gb)",
#                              "Area Of Special Scientific Interest (Ni)",
#                              "Special Protection Area (Birds Directive)",
#                              "Site of Community Importance (Habitats Directive)",
#                              "Nature Reserve")) %>%
#   # british national grid
#   sf::st_transform(crs = gbcrs)
# 
# # template raster - 1km
# # grid reference spatial locations
# rast_temp <- raster::rasterFromXYZ(dplyr::select(gr_ref, x = easting, y = northing), crs = crs(pa))
# 
# # pa coverage at start of time period
# pa_pr_1990 <- pa_pr(pa = pa, rast_temp = rast_temp, gr_ref = gr_ref, yr = 1990)
# 
# 
# # pa coverage at end of time period
# pa_pr_2018 <- pa_pr(pa = pa, rast_temp = rast_temp, gr_ref = gr_ref, yr = 2018)

##### Data #####

## Here we load the protected area data and classify sites as protected or unprotected

# # protected area status
# 
# pa_pr_1990 <- pa_pr_1990 %>%
#   # filter out NAs
#   dplyr::filter(!is.na(grid_ref))
#
# pa_pr_2018 <- pa_pr_2018 %>%
#   # filter out NAs
#   dplyr::filter(!is.na(grid_ref))
# 
# # protected sites
# pro <- pa_pr_1990 %>%
#   # 10 % protection threshold
#   dplyr::filter(pa_prop > 0.1) %>%
#   # protected or not
#   dplyr::mutate(prot = "pa")
# 
# # unprotected sites
# unp <- pa_pr_2018 %>%
#   # 1 % protection threshold
#   dplyr::filter(pa_prop < 0.01) %>%
#   # protected or not
#   dplyr::mutate(prot = "unp")
# 
# # protected status
# pro_stat <- dplyr::bind_rows(pro, unp)
# 
# # save: pro_stat
# saveRDS(pro_stat, "data/pro_stat.rds")

##### Occupancy #####

## First we extract the occupancy data for all taxonomic groups we are interested in

# Taxonomic group(s)
taxa <- c("Ants", "Bees", "Hoverflies", "Ladybirds", "Spiders", "Wasps")
ntax <- length(taxa)

# regions <- c("GB", "pa", "unp")
# nregions <- length(regions)
# 
# vers <- c("2021_Francesca_bwars_rerun", "2021_Francesca_bwars_rerun", "2021_Francesca", "2021_Francesca", "2021_Francesca_marlog_rerun", "2021_Francesca_bwars_rerun")
# 
# # create roster - run for multiple groups & multiple regions
# # ignore warning
# roster <- wrappeR::createRoster(index = 1:(nregions * ntax),
#                        modPath = "/data-s3/occmods/",
#                        metaPath = "/data-s3/metadata/",
#                        ver = rep(vers, each = nregions),
#                        indicator = "all",
#                        region = rep(regions, times = ntax),
#                        nSamps = 999,
#                        minObs = 1,
#                        scaleObs = "region",
#                        write = TRUE,
#                        outPath = "~/proarea/data/filtered/",
#                        group = rep(taxa, each = nregions),
#                        t0 = startyear,
#                        tn = endyear)
# 
# # filtered data
# filt_tax <- lapply(roster, wrappeR::applySamp, parallel = FALSE)

##### Model threshold quality #####

# ## observation metadata
# 
# meta_pa <- lapply(taxa, function(x) {
#   
#   load_rdata(paste0("data/filtered/", x, "_all_pa_samp.rdata")) %>% 
#     .[[2]] %>% 
#     dplyr::mutate(tax_grp = x)
#   
# }) %>% 
#   dplyr::bind_rows(.) 
# 
# pass_keep <- meta_pa %>% 
#   dplyr::filter(rot_EqualWt_r_pa == TRUE) %>% 
#   .$Species_r_pa

##### Occupancy estimates #####

## Here we read in prepared occupancy data

# # read in filtered data
# occ_read <- function(x, reg) {
#   
#   load_rdata(paste0("data/filtered/", x, "_all_", reg, "_samp.rdata")) %>% 
#     .[[1]] %>% 
#     dplyr::mutate(grp = x)
#   
# }
# 
# # focal years 1990 to 2018
# foc_yrs <- paste0("year_", startyear:endyear)
# 
# ## pa
# 
# occ_pa_prep <- lapply(taxa, function(x) occ_read(x = x, reg = "pa")) %>%
#   dplyr::bind_rows(.) %>%
#   # trim to focal years 1990 to 2018
#   dplyr::select(foc_yrs, iteration, species, grp) %>%
#   # add region aggregate name
#   dplyr::mutate(prot = "pa")
# 
# spp_pa <- data.frame(spp = unique(occ_pa_prep$species))
# 
# ## unp
# 
# occ_unp_prep <- lapply(taxa, function(x) occ_read(x = x, reg = "unp")) %>%
#   dplyr::bind_rows(.) %>%
#   # trim to focal years 1990 to 2018
#   dplyr::select(foc_yrs, iteration, species, grp)  %>%
#   # add region aggregate name
#   dplyr::mutate(prot = "unp")
# 
# spp_unp <- data.frame(spp = unique(occ_unp_prep$species))
# 
# # unify protected and unprotected species lists
# 
# all_spp <- dplyr::inner_join(spp_pa, spp_unp) %>% 
#   # remove non-predatory ladybirds
#   dplyr::filter(!spp %in% c("halyzia sedecimguttata", "henosepilachna argus", "psyllobora vigintiduopunctata", "subcoccinella vigintiquattuorpunctata", "tytthaspis sedecimpunctata")) %>% 
#   # remove non-natives
#   dplyr::filter(!spp %in% c("harmonia axyridis", "steatoda nobilis"))
# 
# # update to match species lists
# occ_pa_prep <- dplyr::filter(occ_pa_prep, species %in% all_spp$spp)
# occ_unp_prep <- dplyr::filter(occ_unp_prep, species %in% all_spp$spp)
# 
# # remove species that didn't pass model quality threshold
# occ_pa_prep <- dplyr::filter(occ_pa_prep, species %in% pass_keep)
# occ_unp_prep <- dplyr::filter(occ_unp_prep, species %in% pass_keep)
# 
# GB_spp <- occ_pa_prep %>% 
#   dplyr::distinct(species, .keep_all = TRUE) %>% 
#   dplyr::count(grp) %>% 
#   tibble::add_row(grp = "Overall", n = sum(.$n))
# 
# occ_GB <- lapply(taxa, function(x) occ_read(x = x, reg = "GB")) %>%
#   dplyr::bind_rows(.) %>%
#   # trim to focal years 1990 to 2018
#   dplyr::select(foc_yrs, iteration, species, grp)  %>%
#   # add region aggregate name
#   dplyr::mutate(prot = "GB") %>% 
#   # match species list
#   dplyr::filter(species %in% occ_pa_prep$species) 

#### Functional groups ####

# prepared occupancy data
occ_pa_prep <- readRDS("data/occ_pa_prep.rds")
occ_unp_prep <- readRDS("data/occ_unp_prep.rds")

# duplicate hoverflies as there are primary contributors to two ecosystem functions
hov_pa <- occ_pa_prep %>% 
  dplyr::filter(grp == "Hoverflies") %>% 
  dplyr::mutate(grp = "Hoverflies2")

hov_unp <- occ_unp_prep %>% 
  dplyr::filter(grp == "Hoverflies") %>% 
  dplyr::mutate(grp = "Hoverflies2")

# ecosystem functions
func_lookup <- data.frame(grp = c("Ants", "Hoverflies", "Ladybirds", "Spiders", "Wasps", "Bees", "Hoverflies2"), func = c(rep("Predators", 5), rep("Pollinators", 2)))

occ_pa_func <- occ_pa_prep %>% 
  dplyr::bind_rows(hov_pa) %>% 
  dplyr::left_join(func_lookup, by = "grp") %>% 
  dplyr::rename(tax_grp = grp, grp = func)

occ_unp_func <- occ_unp_prep %>% 
  dplyr::bind_rows(hov_unp) %>% 
  dplyr::left_join(func_lookup, by = "grp") %>% 
  dplyr::rename(tax_grp = grp, grp = func)

occ_pa <- dplyr::bind_rows(occ_pa_func, occ_pa_prep %>% dplyr::mutate(grp = "Overall"))

occ_unp <- dplyr::bind_rows(occ_unp_func, occ_unp_prep %>% dplyr::mutate(grp = "Overall"))

GB_func_spp <- occ_pa %>% 
  dplyr::distinct(species, grp, .keep_all = TRUE) %>% 
  dplyr::count(grp)

# tidy up
rm(hov_pa, hov_unp, occ_pa_func, occ_pa_prep, occ_unp_func, occ_unp_prep)

#### Species richness ####

## Here we calculate species richness across and per year

# function to calculate species richness per year and across years
sprich <- function(occ_df) {
  
  spr_df <- occ_df %>%
    tidyr::gather(year, occ, dplyr::starts_with("year_")) %>%
    dplyr::group_by(year, grp, iteration, prot) %>%
    # species richness per year
    dplyr::summarise(spr_yr = sum(occ, na.rm = TRUE)) %>%
    dplyr::group_by(grp, iteration, prot) %>%
    # mean species richness across years
    dplyr::summarise(spr = mean(spr_yr, na.rm = TRUE)) %>% 
    # tidy protected area status
    dplyr::mutate(prot = dplyr::recode(prot, pa = "Protected", unp = "Unprotected")) %>% 
    dplyr::mutate(prot = factor(prot, levels = c("Unprotected", "Protected")))
  
}

# calculate species richness per year and across years
spr_pa <- sprich(occ_df = occ_pa) %>% 
  dplyr::left_join(GB_func_spp, by = "grp")
spr_unp <- sprich(occ_df = occ_unp) %>% 
  dplyr::left_join(GB_func_spp, by = "grp")

# combine protected and unprotected areas
spr_comb <- dplyr::bind_rows(dplyr::select(spr_pa, grp, iteration, prot, spr), dplyr::select(spr_unp, grp, iteration, prot, spr)) %>% 
  dplyr::ungroup()

# function to calculate difference between protected and unprotected areas
sprich_diff <- function(spr_df_pa, spr_df_unp) {
  
  spr_df <- dplyr::left_join(dplyr::select(spr_df_pa, grp, iteration, spr, pspr), dplyr::select(spr_df_unp, grp, iteration, spr, pspr), by = c("grp", "iteration")) %>% 
    dplyr::mutate(spr_diff = spr.x - spr.y,
                  pspr_diff = pspr.x - pspr.y)
  
  # average across iterations
  spr_av <- spr_df %>% 
    dplyr::group_by(grp) %>% 
    dplyr::summarise_at(vars(spr_diff, pspr_diff), list(
      ~mean(.), 
      ~median(.), 
      low_ci = ~HDInterval::hdi(., credMass = ci)[[1]],
      upp_ci = ~HDInterval::hdi(., credMass = ci)[[2]])) %>% 
    # add significance identifier
    dplyr::mutate(spr_sig = dplyr::case_when(
      spr_diff_low_ci > 0 & spr_diff_upp_ci > 0 ~ "spos",
      spr_diff_low_ci < 0 & spr_diff_upp_ci < 0 ~ "sneg",
      TRUE ~ "ns"
    )) %>% 
    dplyr::mutate(pspr_sig = dplyr::case_when(
      pspr_diff_low_ci > 0 & pspr_diff_upp_ci > 0 ~ "spos",
      pspr_diff_low_ci < 0 & pspr_diff_upp_ci < 0 ~ "sneg",
      TRUE ~ "ns"
    )) %>% 
    dplyr::mutate(spr_sig = factor(spr_sig, levels = c("sneg", "ns", "spos"))) %>% 
    dplyr::mutate(pspr_sig = factor(pspr_sig, levels = c("sneg", "ns", "spos")))
  
  spr_df <- dplyr::left_join(spr_df, dplyr::select(spr_av, grp, spr_sig, pspr_sig), by = c("grp"))
  
  return(list(spr_df, spr_av))
  
}

# calculate difference between protected and unprotected areas
spr_diff_out <- sprich_diff(spr_df_pa = spr_pa, spr_df_unp = spr_unp)

# overall
spr_over <- est_plot_func(grp = "Overall", est_df = spr_comb, est_diff = spr_diff_out, vari = "spr", lim = c(-36, 36), axlab = "Species richness")

# View(spr_over[[2]])
# View(spr_over[[3]])
# View(spr_over[[4]])
# View(spr_over[[6]])

# pollinators
spr_poll <- est_plot_func(grp = "Pollinators", est_df = spr_comb, est_diff = spr_diff_out, vari = "spr", lim = c(-14, 14), axlab = "Species richness")

# predators
spr_pred <- est_plot_func(grp = "Predators", est_df = spr_comb, est_diff = spr_diff_out, vari = "spr", lim = c(-31, 31), axlab = "Species richness")

# # save plot
# cowplot::save_plot("outputs/fig2_sprich_over.png", spr_over[[1]], base_height = 6, base_width = 8, dpi = 300)
