# title: pa_pr_fun
# description: function to rasterize protected areas at 100m and extract proportion protected at 1km
# author(s): Robert Cooke (roboke@ceh.ac.uk)
# date created: 28/10/2020
# date modified: 28/10/2020

pa_pr <- function(pa, rast_temp, gr_ref, yr) {
  
  # protected areas filtered to year (start or end year)
  pa_year <- pa %>% 
    dplyr::filter(STATUS_YR <= yr)
  
  # split template raster
  
  # four parts
  # divide number of columns and rows by 2
  r_spl <- raster::aggregate(rast_temp, c((rast_temp@ncols / 2), (rast_temp@nrows / 2))) 
  r_spl <- as(r_spl, 'SpatialPolygons')
  
  # split raster into four parts
  ra <- lapply(seq_along(r_spl), function(x) raster::crop(rast_temp, r_spl[x]))
  
  # get pa cover for each part
  pa_pr_li <- lapply(1:length(ra), function(x) raster::rasterize(pa_year, ra[[x]], getCover = TRUE))
  
  # merge rasters back together
  pa_pr <- do.call(raster::merge, pa_pr_li)
  
  # extract raster data to dataframe
  pa_pr_df <- as.data.frame(pa_pr, xy = TRUE, centroids = TRUE) %>% 
    # rename layer
    dplyr::rename(pa_prop = layer) %>% 
    # replace nans with na
    dplyr::mutate_at(vars(pa_prop), ~replace(., is.nan(.), NA)) %>% 
    # rename columns
    dplyr::rename(easting = x, northing = y) %>%
    # add grid ref and region data
    dplyr::left_join(gr_ref, by = c("easting", "northing")) %>% 
    # reorder columns
    dplyr::select(grid_ref, region, easting, northing, pa_prop)
  
  # return dataframe
  return(pa_pr_df)
  
}