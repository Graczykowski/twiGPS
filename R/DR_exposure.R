#' Density Ranking exposure
#'
#' @description
#' Using DR function from Yen-Chi Chen's repository density_ranking, only in WGS84 with smoothing parameter
#'
#'
#' @param x data frame with lon lat coordinate columns
#' @param day string in date format compatible with date column in x
#' @param cellsize size of raster cell in meters
#' @param smoothing smoothing parameter
#' @param env_data SparRaster object of envirinmental data
#' @param normalize argument if activity data should be normalized to 0-1 values range
#' @param data_extent TODO
#' @param stats statistics calculated
#' @param act_and_env TODO
#'
#' @return list of SpatRaster result and list of statistics
#'
#'
#'
#' @export


DR_exposure = function(x, day=NULL, cellsize=100, smoothing = 0.0007,
                       normalize = TRUE, env_data=NULL,
                       data_extent = NULL, # TODO extent
                       stats=NULL, act_and_env=FALSE){ # TODO act_and_env

  buff_const = 0.05 #percent of extent as bbox

  if(is.null(day)){
    x_points = x
  }
  else{
    x_points = x  |>  dplyr::filter(wearDate == day)
  }

  # get df with coordinates
  coords = x_points |> dplyr::select("lon", "lat")



  # additional buffer
  lon_range_const = (max(coords[[1]]) - min(coords[[1]])) * buff_const
  lat_range_const = (max(coords[[2]]) - min(coords[[2]])) * buff_const

  # extent of data with buffer
  lon_range = c(min(coords[[1]] - lon_range_const ), max(coords[[1]]) + lon_range_const)
  lat_range = c(min(coords[[2]] - lat_range_const), max(coords[[2]]) + lat_range_const)

  # range distance in WGS84
  dist_lon = geosphere::distm(c(lon_range[1], lat_range[1]), c(lon_range[2], lat_range[1]),
                   fun = geosphere::distHaversine)
  dist_lat = geosphere::distm(c(lon_range[1], lat_range[1]), c(lon_range[1], lat_range[2]),
                   fun = geosphere::distHaversine)

  mean_lonlat_dist = ((dist_lon / cellsize) + (dist_lat / cellsize)) / 2
  # mean cellsize lon and lat


  n_res0 = as.integer(mean_lonlat_dist) # cell number

  # calculate DR data
  DR_data = DR(data=coords, kernel= "Gaussian", h=smoothing, n_res=n_res0,
               xlim=lon_range, ylim=lat_range)


  # creating empty raster to insert KDE values
  len_lon = length(DR_data$x_grid)
  len_lat = length(DR_data$y_grid)

  lon_min = min(DR_data$x_grid)
  lon_max = max(DR_data$x_grid)
  lat_min = min(DR_data$y_grid)
  lat_max = max(DR_data$y_grid)

  spat_dr_rast = terra::rast(nrows=len_lat, ncols=len_lon, xmin=lon_min, xmax=lon_max,
                       ymin=lat_min, ymax=lat_max)

  # insert KDE values to raster
  if (normalize == TRUE){
    dr_values = BBmisc::normalize(DR_data$gr_alpha, method = "range", range = c(0, 1),
                                   margin = 1L, on.constant = "quiet")
  } else {
    dr_values = DR_data$gr_alpha
  }
  terra::values(spat_dr_rast) = dr_values
  spat_dr_rast[spat_dr_rast == 0] = NA # NA vals
  spat_dr_rast = terra::flip(spat_dr_rast, direction='vertical') # flip raster


  if (!is.null(env_data)){ # calculate exposure
    env_data_proj = terra::project(env_data, spat_dr_rast)
    env_data_resamp = terra::resample(env_data_proj, spat_dr_rast)
    rast_env_dr = spat_dr_rast * env_data_resamp

    # calculate env output
    output = output_calc(spat_dr_rast, env_rast = rast_env_dr, stats = stats)
  } else {
    # calculate activity output
    output = output_calc(spat_dr_rast, stats = stats)

  }

  # generate output





  return(output)




}
