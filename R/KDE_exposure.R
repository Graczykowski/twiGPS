#' Kernel Density Estimation exposure (TDA)
#' @description
#' Kernel Density Estimation method activity space and environmental exposure. Using kde() function from TDA package.
#' Bandwidth is a smoothing parameter same as h in TDA::kde(). In order to receive activity space ignore env_data argument.
#'
#' @param x data frame with lon lat coordinates columns
#' @param cellsize positive numeric size of raster cell in meters
#' @param bandwidth positive numeric bandwidth in units of x
#' @param env_data SpatRaster object of environmental data
#' @param normalize boolean argument if activity data should be normalized to 0-1 values range
#' @param data_extent TODO
#' @param start_crs character or terra crs object specifying coordinate reference system of coordinates in x data frame
#' @param end_crs character or terra crs object of coordinate reference system of output
#' @param stats vector of characters statistics to be calculated. See terra::global. "count", "range" and "area" additionally added.
#'
#' @return list of SpatRaster result and data frame of statistics
#'
#' @examples
#'
#' statistics = c("count", "area", "min", "max", "range", "mean", "std", 'sum')
#'
#' data("geolife_sandiego")
#'
#' # activity space
#' KDE_exposure(x = geolife_sandiego, cellsize = 50, bandwidth = 25,
#'  start_crs = "WGS84", end_crs = "EPSG:32611", stats = statistics)
#'
#' #environmental exposure
#' data("landsat_ndvi")
#' ndvi_data = terra::rast(landsat_ndvi)
#'
#' KDE_exposure(x = geolife_sandiego, cellsize = 50, bandwidth = 25,
#'  env_data = ndvi_data, start_crs = "WGS84",
#'  end_crs = "EPSG:32611", stats = statistics)
#'
#'
#'
#' @export
KDE_exposure = function(x, cellsize=NULL, bandwidth = 200, env_data=NULL,
                    normalize = FALSE, data_extent = NULL, # TODO extent
                    start_crs = "WGS84", end_crs=NULL, stats=NULL,
                    act_and_env=FALSE){ # TODO act_and_env

  buff_const = 0.05 #percent of extent as bbox

  x_proj = start_processing(x,  env_data, data_extent, start_crs, end_crs)

  if (!is.null(env_data)){ # change env_data crs beforehand
    env_data_proj = terra::project(env_data, terra::crs(x_proj))
  }



  # get extent
  extent = terra::ext(x_proj)

  # buffer around extent
  x_const_kde = (terra::xmax(extent) - terra::xmin(extent)) * buff_const
  y_const_kde = (terra::ymax(extent) - terra::ymin(extent)) * buff_const

  # new extent
  new_extent = c(terra::xmin(extent) - x_const_kde, terra::xmax(extent) + x_const_kde,
                 terra::ymin(extent) - y_const_kde, terra::ymax(extent) + y_const_kde)
  #TODO small cellsize disproportion


  # implement cellsize > 0 condition
  if (is.numeric(cellsize)) { # cellsize included

    if (terra::linearUnits(x_proj) == 0){ # crs units in degrees
      dist_lon = geosphere::distm(c(new_extent[1], new_extent[3]), c(new_extent[2], new_extent[3]),
                                  fun = geosphere::distHaversine)
      dist_lat = geosphere::distm(c(new_extent[1], new_extent[3]), c(new_extent[1], new_extent[4]),
                                  fun = geosphere::distHaversine)

      x_cells = (dist_lon / cellsize) |> as.integer()

      y_cells = (dist_lat / cellsize) |> as.integer()

      grid_rast = terra::rast(nrows=y_cells, ncols=x_cells, extent = new_extent,
                              crs = terra::crs(x_proj))

    } else {
      grid_rast = terra::rast(crs = terra::crs(x_proj), extent = new_extent,
                              resolution = cellsize)
      ## Probably doesnt matter
      # terra::ext(grid_rast) = new_extent # ext before cellsize to avoid cellsize disproportion
      # terra::res(grid_rast) = cellsize #TODO small cellsize disproportion
      ##
    }

  } else if (!is.null(env_data)){ #if incorrect cellsize and env_data exists

    grid_rast = env_data_proj
    terra::ext(grid_rast) = new_extent # ext before cellsize to avoid cellsize disproportion
    terra::res(grid_rast) = terra::res(env_data_proj)
  }

  # number of cells
  x_cells = terra::ncol(grid_rast)
  y_cells = terra::nrow(grid_rast)




  # coordinate seq
  x_seq = seq(new_extent[1], new_extent[2], length.out = x_cells)
  y_seq = seq(new_extent[3], new_extent[4], length.out = y_cells)

  # coords of every cell
  expand_grid = expand.grid(x_seq,y_seq)

  # point coords
  coords = terra::geom(x_proj)[,3:4]


  # when in degrees smoothing parameter not sqrt until resolving the output of kde
  if (terra::linearUnits(x_proj) == 0){
    kde_data = TDA::kde(coords, Grid = expand_grid, h = bandwidth)
  } else {
    # kde with sqrt bandwidth
    kde_data = TDA::kde(coords, Grid = expand_grid, h = sqrt(bandwidth))
  }

  ### PROBABLY DOESNT MATTER
  # # params for empty raster
  # len_x <- length(x_seq)
  # len_y <- length(y_seq)
  #
  # # x, y limits
  # x_min = min(x_seq)
  # x_max = max(x_seq)
  # y_min = min(y_seq)
  # y_max = max(y_seq)
  #
  # # empty rast for kde
  # kde_rast = terra::rast(nrows=len_y, ncols=len_x, xmin=x_min, xmax=x_max, ymin=y_min,
  #                 ymax=y_max, crs = terra::crs(x_proj))
  ###


  kde_rast = grid_rast
  # insert kde values to raster
  terra::values(kde_rast) = kde_data

  if (normalize == TRUE){ # normalize kde values

    rast_minmax = terra::minmax(kde_rast) # minmax
    # calculate normalization to 0-1 range
    kde_rast = (kde_rast - rast_minmax[1,])/(rast_minmax[2,]-rast_minmax[1,])
  }


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
    output = output_calc(kde_env_rast, stats = stats)
  } else {
    # calculate activity output
    output = output_calc(kde_rast, stats = stats)
  }

  return(output)

}
