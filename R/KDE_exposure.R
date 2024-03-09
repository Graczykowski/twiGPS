#' Kernel Density Estimation in Projected Coordinate Systems
#' @description
#' Code working fast with smoothing parameter as bandwidth still to check
#'
#' @param x data frame with lon lat coordinate columns
#' @param day string in date format compatible with date column in x
#' @param cellsize size of raster cell in meters
#' @param bandwidth bandwidth in unit of x
#' @param env_data SpatRaster object of envirinmental data
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
KDE_exposure = function(x, day=NULL, cellsize=100, bandwidth = 200, env_data=NULL,
                    normalize = FALSE, data_extent = NULL, # TODO extent
                    start_crs = "WGS84", end_crs=NULL, stats=NULL,
                    act_and_env=FALSE){ # TODO act_and_env

  buff_const = 0.05 #percent of extent as bbox

  x_proj = start_processing(x, day, env_data, data_extent, start_crs, end_crs)



  # get extent
  extent = terra::ext(x_proj)

  # buffer around extent
  x_const_kde = (terra::xmax(extent) - terra::xmin(extent)) * buff_const
  y_const_kde = (terra::ymax(extent) - terra::ymin(extent)) * buff_const

  # new extent
  new_extent = c(terra::xmin(extent) - x_const_kde, terra::xmax(extent) + x_const_kde,
                 terra::ymin(extent) - y_const_kde, terra::ymax(extent) + y_const_kde)
  #TODO small cellsize disproportion

  if (!is.null(env_data)){ # change env_data crs beforehand
    env_data_proj = terra::project(env_data, terra::crs(x_proj))
  }



  if (is.numeric(cellsize) & cellsize > 0) { # cellsize included

    if (terra::linearUnits(x_proj) == 0){ # crs units in degrees
      dist_lon = geosphere::distm(c(new_extent[1], new_extent[3]), c(new_extent[2], new_extent[3]),
                                  fun = geosphere::distHaversine)
      dist_lat = geosphere::distm(c(new_extent[1], new_extent[3]), c(new_extent[1], new_extent[4]),
                                  fun = geosphere::distHaversine)

      x_cells = (dist_lon / cellsize) |> as.integer()

      y_cells = (dist_lat / cellsize) |> as.integer()

    } else {
      grid_rast = terra::rast(crs = terra::crs(x_proj))
      terra::ext(grid_rast) = new_extent # ext before cellsize to avoid cellsize disproportion
      terra::res(grid_rast) = cellsize #TODO small cellsize disproportion
    }

  } else if (!is.null(env_data)){ #if incorrect cellsize and env_data exists

    grid_rast = env_data_proj
    terra::ext(grid_rast) = new_extent # ext before cellsize to avoid cellsize disproportion
    terra::res(grid_rast) = terra::res(env_data_proj)
  }

  # number of cells
  if (exists("grid_rast")){
    x_cells = terra::ncol(grid_rast)
    y_cells = terra::nrow(grid_rast)
  }




  # coordinate seq
  x_seq = seq(new_extent[1], new_extent[2], length.out = x_cells)
  y_seq = seq(new_extent[3], new_extent[4], length.out = y_cells)

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
  kde_rast = terra::rast(nrows=len_y, ncols=len_x, xmin=x_min, xmax=x_max, ymin=y_min,
                  ymax=y_max, crs = terra::crs(x_proj))


  if (normalize == TRUE){ # normalize kde values
    kde_vals = BBmisc::normalize(as.vector(kde_data), method = "range", range = c(0, 1),
                                 margin = 1L, on.constant = "quiet")
  } else {
    kde_vals = kde_data
  }

  # insert kde values to raster
  terra::values(kde_rast) = kde_vals
  kde_rast = terra::subst(kde_rast, from = 0, to = NA)
  kde_rast = terra::flip(kde_rast, direction='vertical') # flip raster



  if (!is.null(env_data)){ # calculate exposure
    # only for one point data check
    # env_max = terra::minmax(env_data)[2]
    # env_data_proj[env_data_proj > env_max] = NA
    ##
    env_data_resamp = terra::resample(env_data_proj, kde_rast)
    ##
    # env_data_resamp[env_data_resamp > env_max] = NA
    ##
    kde_env_rast = kde_rast * env_data_resamp


    # calculate env output
    output = output_calc(kde_rast, env_rast = kde_env_rast, stats = stats)
  } else {
    # calculate activity output
    output = output_calc(kde_rast, stats = stats)
  }

  return(output)

}
