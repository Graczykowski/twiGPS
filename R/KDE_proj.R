#' Kernel Density Estimation in Projected Coordinate Systems
#' @description
#' Code working fast with results to check, working well in projected coordinate systems, WGS84 yet to do
#'
#' @param x data frame with lon lat coordinate columns
#' @param day string in date format compatible with date column in x
#' @param cellsize size of raster cell in meters
#' @param bandwidth bandwidth in meters
#' @param env_data SparRaster object of envirinmental data
#' @param normalize argument if activity data should be normalized to 0-1 values range
#' @param data_extent TODO
#' @param start_crs coordinate system of coordinates in x data frame
#' @param end_crs coordinate system of output
#' @param stats statistics calculated
#' @param act_and_env TODO
#'
#' @return list of SpatRaster result and list of statistics
#'
#'
#'
#'
#' @export





KDE_proj = function(x, day=NULL, cellsize=100, bandwidth = 200, env_data=NULL,
                    normalize = FALSE, data_extent = NULL, # TODO extent
                    start_crs = "WGS84", end_crs=NULL, stats=NULL,
                    act_and_env=FALSE){ # TODO act_and_env

  buff_const = 0.05 #percent of extent as bbox

  x_proj = start_processing(x, day, env_data, data_extent, start_crs, end_crs)


  # get extent
  extent = terra::ext(x_proj)

  # buffer around extent
  x_const_kde = (xmax(extent) - xmin(extent)) * buff_const
  y_const_kde = (ymax(extent) - ymin(extent)) * buff_const

  # new extent
  new_extent = c(xmin(extent) - x_const_kde, xmax(extent) + x_const_kde,
                 ymin(extent) - y_const_kde, ymax(extent) + y_const_kde)



  if (is.numeric(cellsize) & cellsize > 0) { # cellsize included

    grid_rast = terra::rast(crs = terra::crs(x_proj))
    terra::ext(grid_rast) = new_extent # ext before cellsize to avoid cellsize disproportion
    terra::res(grid_rast) = cellsize
  } else if (!is.null(env_data)){ #if incorrect cellsize and env_data exists

    grid_rast = env_data
    terra::ext(grid_rast) = new_extent # ext before cellsize to avoid cellsize disproportion
    terra::res(grid_rast) = res(env_data)
  }

  # coordinate seq
  x_seq = seq(new_extent[1], new_extent[2], length.out = terra::dim(grid_rast)[2])
  y_seq = seq(new_extent[3], new_extent[4], length.out = terra::dim(grid_rast)[1])

  # coords of every cell
  expand_grid = expand.grid(x_seq,y_seq)

  # point coords
  coords = terra::geom(x_proj)[,3:4]

  # kde with sqrt bandwidth
  kde_data = TDA::kde(coords, Grid = expand_grid, h = sqrt(bandwidth))

  # params for empty raster
  len_x <- length(x_seq)
  len_y <- length(y_seq)

  # x, y limits
  x_min = min(x_seq)
  x_max = max(x_seq)
  y_min = min(y_seq)
  y_max = max(y_seq)

  # empty rast for kde
  kde_rast = rast(nrows=len_y, ncols=len_x, xmin=x_min, xmax=x_max, ymin=y_min,
                  ymax=y_max, crs = crs(x_proj))


  if (normalize == TRUE){ # normalize kde values
    kde_vals = BBmisc::normalize(kde_rast, method = "range", range = c(0, 1),
                                 margin = 1L, on.constant = "quiet")
  } else {
    kde_vals = kde_data
  }

  # insert kde values to raster
  terra::values(kde_rast) = kde_vals
  kde_rast = terra::subst(kde_rast, from = 0, to = NA)



  if (!is.null(env_data)){ # calculate exposure
    env_data_proj = terra::project(env_data, kde_rast)
    env_data_resamp = terra::resample(env_data_proj, kde_rast)
    kde_env_rast = kde_rast * env_data_resamp


    # calculate env output
    output = output_calc(kde_rast, env_rast = kde_env_rast, stats = stats)
  } else {
    # calculate activity output
    output = output_calc(kde_rast, stats = stats)
  }

  return(output)

}
