#' Kernel Density Estimation in WGS84
#' @description
#' Code similar to DR code but using KDE in WGS84, fast but with smoothing parameter
#'
#' @param x data frame with lon lat coordinate columns
#' @param day string in date format compatible with date column in x
#' @param cellsize size of raster cell in meters
#' @param smoothing smoothing parameter
#' @param env_data SparRaster object of envirinmental data
#' @param data_extent TODO
#' @param stats statistics calculated
#' @param act_and_env TODO
#'
#' @return list of SpatRaster result and list of statistics
#'
#'
#' @export
KDE_wgs = function(x, day=NULL, cellsize=100, smoothing = 0.0007, env_data=NULL,
                   data_extent = NULL, # TODO extent
                   stats=NULL, act_and_env=FALSE){ # TODO act_and_env

  buff_const = 0.05 #percent of extent as bbox

  if(is.null(day)){
    x_points = x
  }
  else{
    x_points = x %>% dplyr::filter(wearDate == day)
  }

  # get df with coordinates
  coords = x_points %>% dplyr::select("lon", "lat")



  # additional buffer
  lon_range_const = (max(coords[[1]]) - min(coords[[1]])) * buff_const
  lat_range_const = (max(coords[[2]]) - min(coords[[2]])) * buff_const

  # extent of data with buffer
  lon_range = c(min(coords[[1]] - lon_range_const ), max(coords[[1]]) + lon_range_const)
  lat_range = c(min(coords[[2]] - lat_range_const), max(coords[[2]]) + lat_range_const)

  # range distance in WGS84
  dist_lon = distm(c(lon_range[1], lat_range[1]), c(lon_range[2], lat_range[1]),
                 fun = geosphere::distHaversine)
  dist_lat = distm(c(lon_range[1], lat_range[1]), c(lon_range[1], lat_range[2]),
                 fun = geosphere::distHaversine)

  mean_lonlat_dist = ((dist_lon / cellsize) + (dist_lat / cellsize)) / 2
  # mean cellsize lon and lat


  n_res0 = as.integer(mean_xy_dist) # cell number

  # creating grid for KDE
  lon_seq = seq(from=lon_range[1], to=lon_range[2], length.out=n_res0)
  lat_seq = seq(from=lat_range[1], to=lat_range[2], length.out=n_res0)
  gr1 = expand.grid(lon_seq,lat_seq)
  names(gr1) = c("x", "y")


  # KDE
  kde_vals = TDA::kde(X = coords, Grid = gr1, h = smoothing)

  # creating empty raster to insert KDE values
  len_lon = length(lon_seq)
  len_lat = length(lat_seq)

  lon_min = min(lon_seq)
  lon_max = max(lon_seq)
  lat_min = min(lat_seq)
  lat_max = max(lat_seq)

  kde_wgs_rast = rast(nrows=len_lat, ncols=len_lon, xmin=lon_min, xmax=lon_max,
                       ymin=lat_min, ymax=lat_max)

  # insert KDE values to raster
  terra::values(kde_wgs_rast) = kde_vals
  kde_wgs_rast[kde_wgs_rast == 0] = NA # NA vals
  kde_wgs_rast = terra::flip(kde_wgs_rast, direction='vertical') # flip raster


  if (!is.null(env_data)){ # calculate exposure
    env_kde_proj = terra::project(env_data, kde_wgs_rast)
    env_kde_resamp = terra::resample(env_kde_proj, kde_wgs_rast)
    rast_env_kde_wgs = kde_wgs_rast * env_kde_resamp


    # calculate env output
    output = output_calc(act_rast = kde_wgs_rast, env_rast = rast_env_kde_wgs, stats = stats)
  } else {
    # calculate activity output
    output = output_calc(act_rast = kde_wgs_rast, stats = stats)
  }



  return(output)



}
