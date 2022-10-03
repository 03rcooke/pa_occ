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
library(maps) # maps: map data
library(mapproj) # mapproj: map manipulation
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

# # install gghalves from github
# remotes::install_github('erocoar/gghalves')
library(gghalves)

# # commit packages to renv library
# renv::snapshot()

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
  
  spr_df <- dplyr::left_join(dplyr::select(spr_df_pa, grp, iteration, spr), dplyr::select(spr_df_unp, grp, iteration, spr), by = c("grp", "iteration")) %>% 
    dplyr::mutate(spr_diff = spr.x - spr.y)
  
  # average across iterations
  spr_av <- spr_df %>% 
    dplyr::group_by(grp) %>% 
    dplyr::summarise_at(vars(spr_diff), list(
      spr_diff_mean = ~mean(.), 
      spr_diff_median = ~median(.), 
      spr_diff_low_ci = ~HDInterval::hdi(., credMass = ci)[[1]],
      spr_diff_upp_ci = ~HDInterval::hdi(., credMass = ci)[[2]]))
  
  spr_df <- dplyr::left_join(spr_df, dplyr::select(spr_av, grp), by = c("grp"))
  
  return(list(spr_df, spr_av))
  
}

# calculate difference between protected and unprotected areas
spr_diff_out <- sprich_diff(spr_df_pa = spr_pa, spr_df_unp = spr_unp)

est_plot_func <- function(grp, est_df, est_diff, vari, lim = NULL, axlab) {
  
  perc_diff_df <- NULL
  
  low_lim <- lim[[1]]
  upp_lim <- lim[[2]]
  
  vari_diff <- paste0(vari, "_diff")
  
  df <- est_df %>% 
    dplyr::filter(grp == !!grp)
  
  diff_df <- est_diff[[1]] %>% 
    dplyr::filter(grp == !!grp) 
  
  if(vari == "spr") {
    
    vari_y <- paste0(vari, ".y")
    
    diff_df <-  diff_df %>% 
      # percentage difference
      dplyr::mutate(perc_diff = (.data[[vari_diff]] / .data[[vari_y]]) * 100)
    
    perc_diff_df <- diff_df %>% 
      dplyr::ungroup() %>% 
      dplyr::summarise_at(vars(perc_diff), list(meds = ~median(.),
                                                low_ci = ~HDInterval::hdi(., credMass = ci)[[1]],
                                                upp_ci = ~HDInterval::hdi(., credMass = ci)[[2]]))
    
  }
  
  diff_av <- est_diff[[2]] %>% 
    dplyr::filter(grp == !!grp)
  
  numb <- GB_func_spp %>% 
    dplyr::filter(grp == !!grp)
  
  # relative effect size
  eff <- effsize::cohen.d(as.formula(paste0(vari, " ~ prot")), data = df, hedges.correction = TRUE)
  
  eff_text <- paste0(sprintf("%.1f", -eff$estimate), " [", sprintf("%.1f", -eff$conf.int[[2]]), ", ", sprintf("%.1f", -eff$conf.int[[1]]), "]")
  
  meds <- df %>% 
    dplyr::group_by(prot) %>% 
    dplyr::summarise_at(vars(vari), list(meds = ~median(.),
                                         low_ci = ~HDInterval::hdi(., credMass = ci)[[1]],
                                         upp_ci = ~HDInterval::hdi(., credMass = ci)[[2]])) %>% 
    dplyr::mutate(diff = c(0, diff_av[ , paste0(vari_diff, "_median")][[1]])) %>% 
    dplyr::mutate(colr = c("#EE7733", "#0077BB"))
  
  est_plot <- ggplot2::ggplot(df, ggplot2::aes(x = prot, y = .data[[vari]])) +
    # median lines
    ggplot2::geom_hline(data = meds, ggplot2::aes(yintercept = meds, colour = prot), lty = 2, lwd = 1) +
    # jittered points
    ggplot2::geom_jitter(ggplot2::aes(colour = prot), alpha = 0.4)  +
    # scales
    ggplot2::scale_y_continuous(name = axlab, limits = c(low_lim + meds[1,2][[1]], upp_lim + meds[1,2][[1]])) +
    ggplot2::scale_x_discrete(labels = c("Unprotected\n ", "Protected\n ")) +
    ggplot2::scale_colour_manual(values = c("#EE7733", "#0077BB")) +
    ggplot2::theme(legend.position = "none",
                   axis.title.x = ggplot2::element_blank())
  
  est_diff_plot <- ggplot2::ggplot(diff_df, ggplot2::aes(x = "", y = .data[[vari_diff]])) +
    # median lines
    ggplot2::geom_hline(data = meds, ggplot2::aes(yintercept = diff), colour = meds$colr, lty = 2, lwd = 1) +
    # half violin
    gghalves::geom_half_violin(fill = "grey50", side = "r", colour = NA, scale = "count", alpha = 0.6) +
    # 95% CrI
    ggplot2::geom_linerange(data = diff_av, ggplot2::aes(ymax = .data[[paste0(vari_diff, "_upp_ci")]], ymin = .data[[paste0(vari_diff, "_low_ci")]], x = ""), colour = "grey50", size = 1, inherit.aes = FALSE) +
    # median
    ggplot2::geom_point(data = diff_av, ggplot2::aes(y = .data[[paste0(vari_diff, "_median")]]), colour = "grey50", size = 10, shape = "-") +
    # relative effect size
    ggplot2::annotate("text", x = 1, y = min(diff_df[vari_diff]), label = eff_text) +
    # scales
    ggplot2::scale_y_continuous(name = "Difference", position = "right", limits = c(low_lim, upp_lim)) +
    ggplot2::scale_x_discrete(labels = c("Protected minus\nunprotected")) +
    ggplot2::theme(legend.position = "none",
                   axis.title.x = ggplot2::element_blank())
  
  est_prep_plot <- cowplot::plot_grid(est_plot, est_diff_plot, rel_widths = c(1, 0.8))
  
  return(list(est_prep_plot,
              meds,
              diff_av,
              eff,
              est_diff_plot,
              perc_diff_df))
  
}

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
# cowplot::save_plot("outputs/fig_2_sprich_over.png", spr_over[[1]], base_height = 6, base_width = 8, dpi = 300)

#### Multi-species trends ####

## Here we calculate geometric mean occupancy

# log of 0 is undefined
nudgeOcc <- function(x, nudgeFac = 0.0001) {
  x[x == 0] <- nudgeFac
  return(x)
}

# function to calculate geometric mean
geomean <- function(x) exp(mean(log(nudgeOcc(x))))

# function to calculate multi-species indicator per iteration
msi <- function(occ_df) {
  
  msi_df <- occ_df %>%
    tidyr::gather(year, occ, dplyr::starts_with("year_")) %>%
    dplyr::group_by(year, grp, iteration, prot) %>%
    # geometric mean
    dplyr::summarise(gm = geomean(occ)) %>% 
    # tidy protected area status
    dplyr::mutate(prot = dplyr::recode(prot, pa = "Protected", unp = "Unprotected")) %>% 
    dplyr::mutate(prot = factor(prot, levels = c("Protected", "Unprotected"))) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(year = as.numeric(gsub("year_", "", year)))
  
}

# calculate multi-species indicator per iteration
msi_pa <- msi(occ_df = occ_pa)
msi_unp <- msi(occ_df = occ_unp)

# function to summarise multi-species indicator across iterations
msi_sum <- function(msi_df) {
  
  # average across iterations
  msi_av <- msi_df %>% 
    dplyr::group_by(prot, grp, year) %>% 
    dplyr::summarise_at(vars(gm), list(
      ~mean(.), 
      ~median(.), 
      low_ci = ~HDInterval::hdi(., credMass = ci)[[1]],
      upp_ci = ~HDInterval::hdi(., credMass = ci)[[2]])) %>%
    dplyr::mutate_at(vars(median), .funs = list(ind = ~ . / (ifelse(!is.na(dplyr::first(.)), dplyr::first(.), dplyr::nth(., 2)) / 100))) %>% 
    dplyr::mutate_at(vars(low_ci, upp_ci), .funs = list(ind = ~ . / (ifelse(!is.na(dplyr::first(median)), dplyr::first(median), dplyr::nth(median, 2)) / 100))) %>% 
    dplyr::ungroup()
  
}

# summarise multi-species indicator across iterations
msi_sum_pa <- msi_sum(msi_df = msi_pa)
msi_sum_unp <- msi_sum(msi_df = msi_unp)

# combine protected and unprotected areas
msi_sum_comb <- dplyr::bind_rows(msi_sum_pa, msi_sum_unp)

# function to plot multi-species trends
trend_plot_func <- function(grp, msi_sum, lim = NULL, leg = FALSE, leg_coord = NULL, axlab = "Geometric mean occupancy") {
  
  low_lim <- lim[[1]]
  upp_lim <- lim[[2]]
  
  numb <- GB_func_spp %>% 
    dplyr::filter(grp == !!grp)
  
  msi <- msi_sum %>% 
    dplyr::filter(grp == !!grp)
  
  refe <- dplyr::filter(msi, year == 1990)
  
  trend_plot <- ggplot(data = msi, aes(x = year, y = median, group = prot)) +
    geom_hline(data = refe, aes(colour = prot, yintercept = median), lty = 2) +
    # 95% CrI
    geom_ribbon(aes(ymin = low_ci, ymax = upp_ci, fill = prot), alpha = 0.6) +
    geom_line(aes(y = median, colour = prot)) +
    geom_point(aes(colour = prot), size = 0.6) +
    # scales
    scale_y_continuous(name = axlab, limits = c(low_lim, upp_lim)) +
    scale_x_continuous(name = "Year", breaks = scales::pretty_breaks(6), expand = c(0,0)) +
    scale_colour_manual(name = "", values = c("#0077BB", "#EE7733")) +
    scale_fill_manual(name = "", values = c("#5DC4FF", "#F6BB99")) +
    theme(legend.position = "none")
  
  if(leg == TRUE) trend_plot <- trend_plot + 
    annotate("text", x = leg_coord[1], y = leg_coord[2], label = "Protected", colour = "#0077BB", size = 6) + 
    annotate("text", x = leg_coord[3], y = leg_coord[4], label = "Unprotected", colour = "#EE7733", size = 6) 
  
  return(trend_plot)
  
}

# overall
msi_over <- trend_plot_func(grp = "Overall", msi_sum = msi_sum_comb, lim = c(0.043, 0.082), leg = TRUE, leg_coord = c(2015, 0.076, 2015, 0.053))

# pollinators
msi_poll <- trend_plot_func(grp = "Pollinators", msi_sum = msi_sum_comb, lim = c(0.037, 0.087), axlab = "Geometric\nmean occupancy")

# predators
msi_pred <- trend_plot_func(grp = "Predators", msi_sum = msi_sum_comb, lim = c(0.044, 0.084), axlab = "Geometric\nmean occupancy")

# Trend change

## Here we calculate annual growth rate between the first and last year

# function to calculate annual growth rate per species
trend_change <- function(occ_df, taxa) {
  
  out <- lapply(taxa, function(x) {
    
    spp_df <- occ_df %>% 
      dplyr::filter(grp == x) %>% 
      # only non-na columns
      dplyr::select(which(colSums(is.na(.)) < 1)) %>% 
      # select first and last year
      dplyr::select(1, (ncol(.) - 4), iteration, species, prot, grp) %>%
      # record first year
      dplyr::mutate(first_year = colnames(.)[1]) %>%
      dplyr::mutate(first_year = as.numeric(gsub("year_", "", first_year))) %>%
      # record last year
      dplyr::mutate(last_year = colnames(.)[2]) %>%
      dplyr::mutate(last_year = as.numeric(gsub("year_", "", last_year))) %>%
      # set column names
      dplyr::rename(occ_first = 1, occ_last = 2) %>% 
      # number of years
      dplyr::mutate(nyr = last_year - first_year) %>%
      # replace zero with very small number - can't divide by zero
      dplyr::mutate(occ_first = ifelse(occ_first == 0, 1e-07, occ_first)) %>%
      # annual growth rate
      dplyr::mutate(change = (((occ_last / occ_first) ^ (1 / nyr)) - 1) * 100)
    
  }) %>% 
    # bind across taxonomic groups
    dplyr::bind_rows(.)
  
}

# calculate annual growth rate per species
trend_pa <- trend_change(occ_df = dplyr::select(occ_pa, -tax_grp), taxa = unique(occ_unp$grp))
trend_unp <- trend_change(occ_df = dplyr::select(occ_unp, -tax_grp), taxa = unique(occ_unp$grp))

# function to summarise annual growth rate per group
tax_change <- function(trend_df) {
  
  chng_tax <- trend_df %>% 
    dplyr::group_by(grp, iteration) %>% 
    dplyr::summarise(change = mean(change))
}

# summarise annual growth rate per group
chng_pa <- tax_change(trend_df = trend_pa)
chng_unp <- tax_change(trend_df = trend_unp)

# combine protected and unprotected areas
chng_comb <- dplyr::bind_rows(
  dplyr::select(chng_pa, grp, iteration, change) %>% 
    dplyr::mutate(prot = "Protected"), 
  dplyr::select(chng_unp, grp, iteration, change) %>% 
    dplyr::mutate(prot = "Unprotected")) %>% 
  dplyr::mutate(prot = factor(prot, levels = c("Unprotected", "Protected")))

# function to calculate difference in growth rates between protected and unprotected areas
change_diff <- function(chng_df_pa, chng_df_unp) {
  
  # combine protected and unprotected
  chng_df <- dplyr::left_join(chng_df_pa, chng_df_unp, by = c("grp", "iteration")) %>% 
    dplyr::rename(change_pa = change.x, change_unp = change.y) %>% 
    dplyr::mutate(change_diff = change_pa - change_unp) %>% 
    dplyr::ungroup()
  
  # average across iterations
  chng_av <- chng_df %>% 
    dplyr::group_by(grp) %>% 
    dplyr::summarise_at(vars(change_diff), list(
      change_diff_mean = ~mean(.), 
      change_diff_median = ~median(.), 
      change_diff_low_ci = ~HDInterval::hdi(., credMass = ci)[[1]],
      change_diff_upp_ci = ~HDInterval::hdi(., credMass = ci)[[2]]))
  
  chng_df <- dplyr::left_join(chng_df, dplyr::select(chng_av, grp), by = c("grp"))
  
  return(list(chng_df, chng_av))
  
}

# calculate difference in growth rates between protected and unprotected areas
chng_diff_out <- change_diff(chng_df_pa = chng_pa, chng_df_unp = chng_unp)

# overall
chng_over <- est_plot_func(grp = "Overall", est_df = chng_comb, est_diff = chng_diff_out, vari = "change", lim = c(-1.8, 1.8), axlab = "Growth rate")

# View(chng_over[[2]])
# View(chng_over[[3]])
# View(chng_over[[4]])

# pollinators
chng_poll <- est_plot_func(grp = "Pollinators", est_df = chng_comb, est_diff = chng_diff_out, vari = "change", lim = c(-2.8, 2.8), axlab = "Growth rate")

# predators
chng_pred <- est_plot_func(grp = "Predators", est_df = chng_comb, est_diff = chng_diff_out, vari = "change", lim = c(-2.4, 2.4), axlab = "Growth rate")

# combine msi and growth rate plots
chng_comb_over <- cowplot::plot_grid(msi_over, chng_over[[1]], nrow = 2, labels = "AUTO")

# # save plot
# cowplot::save_plot("outputs/fig_3_chng_comb_over.png", chng_comb_over, base_height = 11, base_width = 9, dpi = 300)

#### Temporal beta diversity ####

beta_tbi_func <- function(grp, it, occ_df_pa, occ_df_unp) {

  # select useful data
  df_pa <- dplyr::select(occ_df_pa, dplyr::starts_with("year"), iteration, species, grp)

  df_unp <- dplyr::select(occ_df_unp, dplyr::starts_with("year"), iteration, species, grp)

  # prepare data as wide species dataframe for first and last year
  bprep_pa <- df_pa %>%
    dplyr::filter(grp == !!grp, iteration == it) %>%
    # only non-na columns
    dplyr::select(which(colSums(is.na(.)) < 1)) %>%
    # select first and last year
    dplyr::select(1, (ncol(.) - 4), species) %>%
    # make truly long
    tidyr::gather(year, occ, dplyr::starts_with("year")) %>%
    # spread to wide
    tidyr::spread(species, occ) %>%
    dplyr::select(-year)

  bprep_unp <- df_unp %>%
    dplyr::filter(grp == !!grp, iteration == it) %>%
    # only non-na columns
    dplyr::select(which(colSums(is.na(.)) < 1)) %>%
    # select first and last year
    dplyr::select(1, (ncol(.) - 4), species) %>%
    # make truly long
    tidyr::gather(year, occ, dplyr::starts_with("year")) %>%
    # spread to wide
    tidyr::spread(species, occ) %>%
    dplyr::select(-year)

  # time 1
  t1 <- bind_rows(bprep_pa[1,], bprep_unp[1,])

  # time 2
  t2 <- bind_rows(bprep_pa[2,], bprep_unp[2,])

  # calculate tbi
  tbi <- adespatial::TBI(t1, t2)

  # compile data in easy format
  tbi_df <- data.frame(tbi_val = tbi$TBI, loss = tbi$BCD.mat$`B/(2A+B+C)`, gain = tbi$BCD.mat$`C/(2A+B+C)`, dir = tbi$BCD.mat$Change, prot = c("Protected", "Unprotected"), iteration = it, grp = grp)

  return(tbi_df)

}

# # calculate temporal beta diversity per iteration
# # takes a while so save output
# beta_tbi <- lapply(unique(occ_pa$grp), function(xgrp) {
# 
#   out <- lapply(1:999, function(x) beta_tbi_func(grp = xgrp, it = x, occ_df_pa = occ_pa, occ_df_unp = occ_unp) %>%
#     # bind across iterations
#   dplyr::bind_rows()
# 
# }) %>%
#   # bind across groups
#   dplyr::bind_rows()
# 
# saveRDS(beta_tbi, "data/beta_tbi.rds")

# beta diversity estimates
beta_tbi <- readRDS("data/beta_tbi.rds") %>% 
  dplyr::mutate(prot = factor(prot, levels = c("Unprotected", "Protected"))) %>% 
  # convert to percentages
  dplyr::mutate_at(vars(tbi_val, loss, gain), ~. * 100)

# contrast beta diversity estimates
beta_contr <- dplyr::left_join(dplyr::filter(beta_tbi, prot == "Protected"), dplyr::filter(beta_tbi, prot == "Unprotected"), by = c("iteration", "grp"))

# function to calculate difference in beta diversity between protected and unprotected areas
beta_diff <- function(beta_df) {
  
  beta_df <- beta_df %>% 
    dplyr::mutate(tbi_val_diff = tbi_val.x - tbi_val.y,
                  loss_diff = loss.x - loss.y,
                  gain_diff = gain.x - gain.y)
  
  # average across iterations
  beta_av <- beta_df %>% 
    dplyr::group_by(grp) %>% 
    dplyr::summarise_at(vars(dplyr::ends_with("diff")), list(
      ~mean(.), 
      ~median(.), 
      low_ci = ~HDInterval::hdi(., credMass = ci)[[1]],
      upp_ci = ~HDInterval::hdi(., credMass = ci)[[2]]))
  
  beta_df <- dplyr::left_join(beta_df, dplyr::select(beta_av, grp), by = c("grp"))
  
  return(list(beta_df, beta_av))
}

# calculate difference in beta diversity between protected and unprotected areas
beta_diff_out <- beta_diff(beta_df = beta_contr)

## total temporal beta diversity

# overall
beta_over <- est_plot_func(grp = "Overall", est_df = beta_tbi, est_diff = beta_diff_out, vari = "tbi_val", lim = c(-2.2, 2.2), axlab = "Temporal beta diversity (%)")

# View(beta_over[[2]])
# View(beta_over[[3]])
# View(beta_over[[4]])

# net change
# overall
nc_over <- beta_tbi %>% 
  dplyr::filter(grp == "Overall") %>% 
  dplyr::mutate(nc = gain - loss) %>% 
  dplyr::group_by(prot) %>% 
  dplyr::summarise_at(vars(nc), list(meds = ~median(.),
                                     low_ci = ~HDInterval::hdi(., credMass = ci)[[1]],
                                     upp_ci = ~HDInterval::hdi(., credMass = ci)[[2]]))

# pollinators
beta_poll <- est_plot_func(grp = "Pollinators", est_df = beta_tbi, est_diff = beta_diff_out, vari = "tbi_val", lim = c(-3.1, 3.1), axlab = "Temporal\nbeta diversity (%)")

# predators
beta_pred <- est_plot_func(grp = "Predators", est_df = beta_tbi, est_diff = beta_diff_out, vari = "tbi_val", lim = c(-2.5, 2.5), axlab = "Temporal\nbeta diversity (%)")

## loss

# overall
loss_over <- est_plot_func(grp = "Overall", est_df = beta_tbi, est_diff = beta_diff_out, vari = "loss", lim = c(-2.7, 2.7), axlab = "Loss (%)")

# View(loss_over[[2]])
# View(loss_over[[3]])
# View(loss_over[[4]])

# pollinators
loss_poll <- est_plot_func(grp = "Pollinators", est_df = beta_tbi, est_diff = beta_diff_out, vari = "loss", lim = c(-2.3, 2.3), axlab = "Loss (%)")

# predators
loss_pred <- est_plot_func(grp = "Predators", est_df = beta_tbi, est_diff = beta_diff_out, vari = "loss", lim = c(-3.8, 3.8), axlab = "Loss (%)")

## gain

# overall
gain_over <- est_plot_func(grp = "Overall", est_df = beta_tbi, est_diff = beta_diff_out, vari = "gain", lim = c(-2.4, 2.4), axlab = "Gain (%)")

# View(gain_over[[2]])
# View(gain_over[[3]])
# View(gain_over[[4]])

# pollinators
gain_poll <- est_plot_func(grp = "Pollinators", est_df = beta_tbi, est_diff = beta_diff_out, vari = "gain", lim = c(-3.1, 3.1), axlab = "Gain (%)")

# predators
gain_pred <- est_plot_func(grp = "Predators", est_df = beta_tbi, est_diff = beta_diff_out, vari = "gain", lim = c(-3.1, 3.1), axlab = "Gain (%)")

# beta diversity plot for overall
beta_over_plot <- cowplot::plot_grid(NULL, beta_over[[1]], NULL, nrow = 1, rel_widths = c(0.4, 1, 0.4), labels = c("", "A", "")) %>% cowplot::plot_grid(., cowplot::plot_grid(loss_over[[1]], gain_over[[1]], nrow = 1, labels = c("B", "C")), nrow = 2)

# # save plot
# cowplot::save_plot("outputs/fig_4_beta_over.png", beta_over_plot, base_height = 8, base_width = 12, dpi = 300)

#### Pollinators combined plots #####

# net change
# pollinators
nc_poll <- beta_tbi %>% 
  dplyr::filter(grp == "Pollinators") %>% 
  dplyr::mutate(nc = gain - loss) %>% 
  dplyr::group_by(prot) %>% 
  dplyr::summarise_at(vars(nc), list(meds = ~median(.),
                                     low_ci = ~HDInterval::hdi(., credMass = ci)[[1]],
                                     upp_ci = ~HDInterval::hdi(., credMass = ci)[[2]]))

# beta diversity
beta_poll_plot <- cowplot::plot_grid(NULL, beta_poll[[1]], NULL, nrow = 1, rel_widths = c(0.4, 1, 0.4), labels = c("", "D", "")) %>% cowplot::plot_grid(., cowplot::plot_grid(loss_poll[[1]], gain_poll[[1]], nrow = 1, labels = c("E", "F")), nrow = 2)

# species richness, trends, beta diversity
poll_plot <- cowplot::plot_grid(cowplot::plot_grid(NULL, spr_poll[[1]], NULL, labels = c("", "A", ""), rel_widths = c(0.4, 1, 0.4), ncol = 3), cowplot::plot_grid(msi_poll, chng_poll[[1]], labels = c("B", "C"), ncol = 2), beta_poll_plot, nrow = 3, rel_heights = c(1, 1, 2))

# # save plot
# cowplot::save_plot("outputs/fig_5_poll.png", poll_plot, base_height = 12, base_width = 10, dpi = 300)

#### Predators combined plots #####

# beta diversity
beta_pred_plot <- cowplot::plot_grid(NULL, beta_pred[[1]], NULL, nrow = 1, rel_widths = c(0.4, 1, 0.4), labels = c("", "D", "")) %>% cowplot::plot_grid(., cowplot::plot_grid(loss_pred[[1]], gain_pred[[1]], nrow = 1, labels = c("E", "F")), nrow = 2)

# species richness, trends, beta diversity
pred_plot <- cowplot::plot_grid(cowplot::plot_grid(NULL, spr_pred[[1]], NULL, labels = c("", "A", ""), rel_widths = c(0.4, 1, 0.4), ncol = 3), cowplot::plot_grid(msi_pred, chng_pred[[1]], labels = c("B", "C"), ncol = 2), beta_pred_plot, nrow = 3, rel_heights = c(1, 1, 2))

# # save plot
# cowplot::save_plot("outputs/appendix_C_fig_s2_pred.png", pred_plot, base_height = 12, base_width = 10, dpi = 300)

#### GB plot ####

# fig 1
# plot of protected vs unprotected grid cells

# sites used
sites_gr <- readRDS("data/sites_gr.rds")

pro_plot <- ggplot2::ggplot(data = sites_gr, aes(x = EASTING, y = NORTHING, fill = prot)) +
  ggplot2::geom_raster() +
  ggplot2::coord_equal() +
  ggplot2::scale_fill_manual(breaks = c("pa", "unp"), labels = c("Protected", "Unprotected"), values = c("#0077BB", "#EE7733", "grey")) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.title = element_blank(),
                 legend.position = "bottom")

euro <- ggplot2::map_data("world") %>% 
  dplyr::mutate(GB = as.factor(ifelse(subregion %in% c("Great Britain", "Wales", "Isle of Wight", "Scotland"), 1, 0)))

euro_plot <- ggplot2::ggplot(euro, aes(x = long, y = lat)) +
  ggplot2::geom_polygon(aes(group = group, fill = GB, colour = GB), size = 0.3) +
  ggplot2::coord_map(project = "orthographic", xlim = c(-12,44), ylim = c(35,70)) +
  ggplot2::scale_fill_manual(values = c("grey", "black")) +
  ggplot2::scale_colour_manual(values = c("darkgrey", "black")) +
  ggplot2::scale_x_continuous(breaks = seq(-140, 60, by = 20))+
  ggplot2::scale_y_continuous(breaks = seq(10, 90, by = 20)) +
  ggplot2::theme_bw() +
  ggplot2::theme(panel.grid.major = element_line(colour = "lightgrey"),
                 panel.border = element_rect(colour = "black"),
                 axis.title = element_blank(),
                 axis.ticks = element_blank(),
                 axis.text = element_blank(),
                 legend.position = "none")

pro_sub_plot <- cowplot::ggdraw(pro_plot) +
  cowplot::draw_plot(euro_plot, 0.63, 0.55, 0.35, 0.35)

# # save: pro_sub_plot
# cowplot::save_plot("outputs/fig_1_pro_sub.png", pro_sub_plot, base_height = 6, base_width = 6, dpi = 300)

#### Number of records ####

# meta_pa <- lapply(taxa, function(x) {
# 
#   load_rdata(paste0("data/filtered/", x, "_all_pa_samp.rdata")) %>%
#     .[[2]] %>%
#     dplyr::mutate(tax_grp = x)
# 
# }) %>%
#   dplyr::bind_rows(.) %>%
#   dplyr::filter(Species_r_pa %in% occ_pa$species)
# 
# rec_pa <- meta_pa %>% 
#   dplyr::group_by(tax_grp) %>% 
#   dplyr::summarise(recs = sum(n_obs_regional_r_pa))
# 
# meta_unp <- lapply(taxa, function(x) {
# 
#   load_rdata(paste0("data/filtered/", x, "_all_unp_samp.rdata")) %>%
#     .[[2]] %>%
#     dplyr::mutate(tax_grp = x)
# 
# }) %>%
#   dplyr::bind_rows(.) %>%
#   dplyr::filter(Species_r_unp %in% occ_unp$species)
# 
# rec_unp <- meta_unp %>%
#   dplyr::group_by(tax_grp) %>%
#   dplyr::summarise(recs = sum(n_obs_regional_r_unp))
# 
# rec_all <- dplyr::bind_cols(rec_pa, rec_unp) %>% 
#   dplyr::mutate(rec_all = recs + recs1)
# 
# saveRDS(rec_all, "data/rec_all.rds")

rec_all <- readRDS("data/rec_all.rds")

sum(rec_all$rec_all)

#### Preprocessing input data for ROBITT - occAssess ####

# # function to convert first letter to uppercase
# firstup <- function(x) {
#   substr(x, 1, 1) <- toupper(substr(x, 1, 1))
#   x
# }
# 
# # function to tidy observations used in models
# # occ_pa and sites_gr needed
# obs <- function(tax_path, tax_grp, filetype, chained = FALSE, min_iter) {
#   
#   spp <- dplyr::filter(occ_pa, tax_grp == !!tax_grp)$species %>% 
#     unique(.) %>% 
#     firstup(.)
#   
#   all_out <- lapply(spp, function(sp) {
#     
#     if(isTRUE(chained)) {
#       
#       out_meta <- load_rdata(paste0(tax_path, sp, "_", min_iter, "_1.rdata"))
#       
#     } else {
#       
#       if(filetype == "rds") {
#         
#         out_meta <- readRDS(paste0(tax_path, sp, ".rds"))
#         
#         }
#     
#       else if(filetype == "rdata") {
#       
#         out_meta <- load_rdata(paste0(tax_path, sp, ".rdata"))
#       
#       }
#     }
#     
#     dat <- out_meta$model$data()
#     
#     dat_df <- data.frame(year = dat$Year, # year
#                          rec = dat$y, # records
#                          site = dat$Site) # sites
# 
#     # site number to grid ref look-up
#     dat_site <- dat_df %>% 
#       dplyr::distinct(site) %>% 
#       dplyr::mutate(grid = unique(out_meta$sites_included[[1]]))
#     
#     region_site_pa <- dat[[paste0("r_", "pa")]][dat$Site]
#     region_site_unp <- dat[[paste0("r_", "unp")]][dat$Site]
#     
#     # sites included within region
#     dat_reg <- cbind(dat_df, region_site_pa, region_site_unp) %>% 
#       dplyr::left_join(dat_site)
#           
#     # subset to temporal window
#     dat_reg_pa <- dat_reg[dat_reg$region_site_pa == 1 & dat_reg$year >= (startyear - (out_meta$min_year - 1)) & dat_reg$year <= (endyear - (out_meta$min_year - 1)), ] %>% 
#       dplyr::filter(rec == 1) %>% 
#       dplyr::mutate(prot = "pa")
#     
#     dat_reg_unp <- dat_reg[dat_reg$region_site_unp == 1 & dat_reg$year >= (startyear - (out_meta$min_year - 1)) & dat_reg$year <= (endyear - (out_meta$min_year - 1)), ] %>% 
#       dplyr::filter(rec == 1) %>% 
#       dplyr::mutate(prot = "unp")
#     
#     dat_comb <- dplyr::bind_rows(dat_reg_pa, dat_reg_unp) %>% 
#       dplyr::mutate(year = (year + 1970) - 1) %>% 
#       dplyr::left_join(dplyr::select(sites_gr, -prot), by = c("grid" = "grid_ref")) %>% 
#       dplyr::mutate(species = tolower(sp)) %>% 
#       dplyr::mutate(tax_grp = tax_grp) %>% 
#       dplyr::select(species, tax_grp, prot, grid, EASTING, NORTHING, year)
#     
#   }) %>% 
#     dplyr::bind_rows()
#   
#   return(all_out)
#   
# }
# 
# ants_obs <- obs(tax_path = "/data-s3/occmods/Ants/occmod_outputs/2021_Francesca_bwars_rerun/", tax_grp = "Ants", filetype = "rdata")
# 
# # saveRDS(ants_obs, "data/ants_obs.rds")
# 
# bees_obs <- obs(tax_path = "/data-s3/occmods/Bees/occmod_outputs/2021_Francesca_bwars_rerun/", tax_grp = "Bees", filetype = "rdata")
# 
# # saveRDS(bees_obs, "data/bees_obs.rds")
# 
# hover_obs <- obs(tax_path = "/data-s3/occmods/Hoverflies/occmod_outputs/2021_Francesca/", tax_grp = "Hoverflies", filetype = "rdata", chained = TRUE, min_iter = 2000)
# 
# # saveRDS(hover_obs, "data/hover_obs.rds")
# 
# lady_obs <- obs(tax_path = "/data-s3/occmods/Ladybirds/occmod_outputs/2021_Francesca/", tax_grp = "Ladybirds", filetype = "rdata", chained = TRUE, min_iter = 8000)
# 
# # saveRDS(lady_obs, "data/lady_obs.rds")
# 
# spider_obs <- obs(tax_path = "/data-s3/occmods/Spiders/occmod_outputs/2021_Francesca_marlog_rerun/", tax_grp = "Spiders", filetype = "rdata")
# 
# # saveRDS(spider_obs, "data/spider_obs.rds")
# 
# wasps_obs <- obs(tax_path = "/data-s3/occmods/Wasps/occmod_outputs/2021_Francesca_bwars_rerun/", tax_grp = "Wasps", filetype = "rdata")
# 
# # saveRDS(wasps_obs, "data/wasps_obs.rds")
# 
# all_obs <- dplyr::bind_rows(ants_obs, bees_obs, hover_obs, lady_obs, spider_obs, wasps_obs)
# 
# # saveRDS(all_obs, "data/all_obs.rds")

all_obs <- readRDS("data/all_obs.rds")

#### rhat ####

# rhat <- function(tax_path, tax_grp, filetype, chained = FALSE, max_iter) {
#   
#   spp <- dplyr::filter(occ_pa, tax_grp == !!tax_grp)$species %>% 
#     unique(.) %>% 
#     firstup(.)
#   
#   bugs <- lapply(spp, function(sp) {
#     
#     if(isTRUE(chained)) {
#       
#       out_dat <- load_rdata(paste0(tax_path, sp, "_", max_iter, "_1.rdata"))
#       
#       if(!is.null(out_dat) & is.null(out_dat$model)) out_dat <- out_dat$out
#       
#     } else {
#       
#       if(filetype == "rds") {
#         
#         out_dat <- readRDS(paste0(tax_path, sp, ".rds"))
#         
#         }
#     
#       else if(filetype == "rdata") {
#       
#         out_dat <- load_rdata(paste0(tax_path, sp, ".rdata"))
#       
#       }
#     }
#     
#     bugs <- out_dat$BUGSoutput$summary %>% 
#       data.frame() %>% 
#       tibble::rownames_to_column("para") %>% 
#       dplyr::filter(stringr::str_detect(para, "psi.fs")) %>% 
#       dplyr::filter(!stringr::str_detect(para, "psi.fs.r_ENGLAND|psi.fs.r_GB|psi.fs.r_NORTHERN_IRELAND|psi.fs.r_SCOTLAND|psi.fs.r_UK|psi.fs.r_WALES|psi.fs.r_high|psi.fs.r_low|psi.fs.r_no_agri|psi.fs.r_pa|psi.fs.r_unp")) %>% 
#       dplyr::mutate(species = tolower(sp)) %>% 
#       dplyr::mutate(tax_grp = tax_grp)
#     
#   }) %>% 
#     dplyr::bind_rows()
#   
#   return(bugs)
#   
# }
# 
# ants_rhat <- rhat(tax_path = "/data-s3/occmods/Ants/occmod_outputs/2021_Francesca_bwars_rerun/", tax_grp = "Ants", filetype = "rdata")
# 
# # saveRDS(ants_rhat, "data/ants_rhat.rds")
# 
# bees_rhat <- rhat(tax_path = "/data-s3/occmods/Bees/occmod_outputs/2021_Francesca_bwars_rerun/", tax_grp = "Bees", filetype = "rdata")
# 
# # saveRDS(bees_rhat, "data/bees_rhat.rds")
# 
# hover_rhat <- rhat(tax_path = "/data-s3/occmods/Hoverflies/occmod_outputs/2021_Francesca/", tax_grp = "Hoverflies", filetype = "rdata", chained = TRUE, max_iter = 32000)
# 
# # saveRDS(hover_rhat, "data/hover_rhat.rds")
# 
# lady_rhat <- rhat(tax_path = "/data-s3/occmods/Ladybirds/occmod_outputs/2021_Francesca/", tax_grp = "Ladybirds", filetype = "rdata", chained = TRUE, max_iter = 32000)
# 
# # saveRDS(lady_rhat, "data/lady_rhat.rds")
# 
# spider_rhat <- rhat(tax_path = "/data-s3/occmods/Spiders/occmod_outputs/2021_Francesca_marlog_rerun/", tax_grp = "Spiders", filetype = "rdata")
# 
# # saveRDS(spider_rhat, "data/spider_rhat.rds")
# 
# wasps_rhat <- rhat(tax_path = "/data-s3/occmods/Wasps/occmod_outputs/2021_Francesca_bwars_rerun/", tax_grp = "Wasps", filetype = "rdata")
# 
# # saveRDS(wasps_rhat, "data/wasps_rhat.rds")
# 
# all_rhat <- dplyr::bind_rows(ants_rhat, bees_rhat, hover_rhat, lady_rhat, spider_rhat, wasps_rhat)
# 
# # saveRDS(all_rhat, "data/all_rhat.rds")

all_rhat <- readRDS("data/all_rhat.rds")

# first and last rhats
fl_rhat <- all_rhat %>% 
  dplyr::group_by(species) %>% 
  slice(c(1,n()))

conv_spp <- fl_rhat %>% 
  dplyr::group_by(species) %>%
  dplyr::filter(all(Rhat < 1.1))

#### Species whose models converged ####

## Rhat < 1.1 for first and last years

# filter to species whose models converged
occ_pa_conv <- dplyr::filter(occ_pa, species %in% conv_spp$species)
occ_unp_conv <- dplyr::filter(occ_unp, species %in% conv_spp$species)

# species richness
spr_pa_conv <- sprich(occ_df = occ_pa_conv) %>% 
  dplyr::left_join(GB_func_spp, by = "grp")
spr_unp_conv <- sprich(occ_df = occ_unp_conv) %>% 
  dplyr::left_join(GB_func_spp, by = "grp")

spr_comb_conv <- dplyr::bind_rows(dplyr::select(spr_pa_conv, grp, iteration, prot, spr), dplyr::select(spr_unp_conv, grp, iteration, prot, spr))

spr_diff_out_conv <- sprich_diff(spr_df_pa = spr_pa_conv, spr_df_unp = spr_unp_conv)

spr_over_conv <- est_plot_func(grp = "Overall", est_df = spr_comb_conv, est_diff = spr_diff_out_conv, vari = "spr", lim = c(-14, 14), axlab = "Species richness")

# multi-species indicator
msi_pa_conv <- msi(occ_df = occ_pa_conv)
msi_unp_conv <- msi(occ_df = occ_unp_conv)

msi_sum_pa_conv <- msi_sum(msi_df = msi_pa_conv)
msi_sum_unp_conv <- msi_sum(msi_df = msi_unp_conv)

msi_sum_comb_conv <- dplyr::bind_rows(msi_sum_pa_conv, msi_sum_unp_conv)

msi_over_conv <- trend_plot_func(grp = "Overall", msi_sum = msi_sum_comb_conv, lim = c(0.055, 0.12), leg = FALSE, axlab = "Geometric\nmean occupancy")

# growth rates
trend_pa_conv <- trend_change(occ_df = dplyr::select(occ_pa_conv, -tax_grp), taxa = unique(occ_unp_conv$grp))
trend_unp_conv <- trend_change(occ_df = dplyr::select(occ_unp_conv, -tax_grp), taxa = unique(occ_unp_conv$grp))

chng_pa_conv <- tax_change(trend_df = trend_pa_conv)
chng_unp_conv <- tax_change(trend_df = trend_unp_conv)

chng_comb_conv <- dplyr::bind_rows(
  dplyr::select(chng_pa_conv, grp, iteration, change) %>% 
    dplyr::mutate(prot = "Protected"), 
  dplyr::select(chng_unp_conv, grp, iteration, change) %>% 
    dplyr::mutate(prot = "Unprotected")) %>% 
  dplyr::mutate(prot = factor(prot, levels = c("Unprotected", "Protected")))

chng_diff_out_conv <- change_diff(chng_df_pa = chng_pa_conv, chng_df_unp = chng_unp_conv)

chng_over_conv <- est_plot_func(grp = "Overall", est_df = chng_comb_conv, est_diff = chng_diff_out_conv, vari = "change", lim = c(-2.8, 2.8), axlab = "Growth rate")

chng_comb_over_conv <- cowplot::plot_grid(msi_over_conv, chng_over_conv[[1]], nrow = 2, labels = "AUTO")

# # beta diversity
# beta_tbi_conv <- lapply(unique(occ_pa_conv$grp), function(xgrp) {
# 
#   out <- lapply(1:999, function(x) beta_tbi_func(grp = xgrp, it = x, occ_df_pa = occ_pa_conv, occ_df_unp = occ_unp_conv)) %>%
#     # bind across iterations
#   dplyr::bind_rows()
# 
# }) %>%
#   # bind across groups
#   dplyr::bind_rows()
# 
# saveRDS(beta_tbi_conv, "data/beta_tbi_conv.rds")

beta_tbi_conv <- readRDS("data/beta_tbi_conv.rds") %>% 
  dplyr::mutate(prot = factor(prot, levels = c("Unprotected", "Protected"))) %>% 
  # convert to percentages
  dplyr::mutate_at(vars(tbi_val, loss, gain), ~. * 100)

beta_contr_conv <- dplyr::left_join(dplyr::filter(beta_tbi_conv, prot == "Protected"), dplyr::filter(beta_tbi_conv, prot == "Unprotected"), by = c("iteration", "grp"))

beta_diff_out_conv <- beta_diff(beta_df = beta_contr_conv)

## total temporal beta diversity

beta_over_conv <- est_plot_func(grp = "Overall", est_df = beta_tbi_conv, est_diff = beta_diff_out_conv, vari = "tbi_val", lim = c(-2.8, 2.8), axlab = "Temporal\nbeta diversity (%)")

## loss

loss_over_conv <- est_plot_func(grp = "Overall", est_df = beta_tbi_conv, est_diff = beta_diff_out_conv, vari = "loss", lim = c(-2.3, 2.3), axlab = "Loss (%)")

# View(loss_over_conv[[2]])
# View(loss_over_conv[[3]])
# View(loss_over_conv[[4]])

## gain

gain_over_conv <- est_plot_func(grp = "Overall", est_df = beta_tbi_conv, est_diff = beta_diff_out_conv, vari = "gain", lim = c(-1.4, 1.4), axlab = "Gain (%)")

# View(gain_over_conv[[2]])
# View(gain_over_conv[[3]])
# View(gain_over_conv[[4]])

beta_over_plot_conv <- cowplot::plot_grid(NULL, beta_over_conv[[1]], NULL, nrow = 1, rel_widths = c(0.4, 1, 0.4), labels = c("", "D", "")) %>% cowplot::plot_grid(., cowplot::plot_grid(loss_over_conv[[1]], gain_over_conv[[1]], nrow = 1, labels = c("E", "F")), nrow = 2)

conv_plot <- cowplot::plot_grid(cowplot::plot_grid(NULL, spr_over_conv[[1]], NULL, labels = c("", "A", ""), rel_widths = c(0.4, 1, 0.4), ncol = 3), cowplot::plot_grid(msi_over_conv, chng_over_conv[[1]], labels = c("B", "C"), ncol = 2), beta_over_plot_conv, nrow = 3, rel_heights = c(1, 1, 2))

# # save plot
# cowplot::save_plot("outputs/appendix_C_fig_s1_conv.png", conv_plot, base_height = 12, base_width = 10, dpi = 300)

