#' Density Ranking exposure in cells degrees unit
#'
#' @description
#' Using DR function from Yen-Chi Chen's repository density_ranking, only in WGS84 with smoothing parameter
#'
#'
#' @param x data frame with lon lat coordinate columns
#' @param day string in date format compatible with date column in x
#' @param cellsize size of raster cell in degrees
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


DR_cellwgs = function(x, day=NULL, cellsize=100, smoothing = 0.0007,
                       normalize = TRUE, env_data=NULL,
                       data_extent = NULL, # TODO extent
                       stats=NULL, act_and_env=FALSE){ # TODO act_and_env

  buff_const = 0.05 #percent of extent as bbox

  x_proj = start_processing(x, day, env_data, data_extent, "WGS84", "WGS84")


  # get extent
  extent = terra::ext(x_proj)

  # buffer around extent
  x_const_dr = (terra::xmax(extent) - terra::xmin(extent)) * buff_const
  y_const_dr = (terra::ymax(extent) - terra::ymin(extent)) * buff_const

  # new extent
  new_extent = c(terra::xmin(extent) - x_const_dr, terra::xmax(extent) + x_const_dr,
                 terra::ymin(extent) - y_const_dr, terra::ymax(extent) + y_const_dr)



  if (is.numeric(cellsize) & cellsize > 0) { # cellsize included

    grid_rast = terra::rast(crs = terra::crs(x_proj))
    terra::ext(grid_rast) = new_extent # ext before cellsize to avoid cellsize disproportion
    terra::res(grid_rast) = cellsize
  }

  # point coords
  coords = terra::geom(x_proj)[,3:4]

  # number of cells
  x_cells = terra::ncol(grid_rast)
  y_cells = terra::nrow(grid_rast)



  # calculate DR data
  dr_data = DR_simple(coords, kernel= "Gaussian", h = smoothing,
                      x_res=x_cells, y_res=y_cells, xlim=new_extent[1:2],
                      ylim=new_extent[3:4])


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
    # only for one point data check
    env_max = terra::minmax(env_data)[2]
    env_data_proj[env_data_proj > env_max] = NA
    ##
    env_data_resamp = terra::resample(env_data_proj, spat_dr_rast)
    ##
    env_data_resamp[env_data_resamp > env_max] = NA
    ##
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
