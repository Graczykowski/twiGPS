#' Kernel Density Estimation exposure (TDA)
#' @description
#' Kernel Density Estimation method activity space and environmental exposure. Using kde() function from TDA package.
#' Bandwidth is a smoothing parameter same as h in TDA::kde(). In order to receive activity space ignore env_data argument.
#'
#' @param x data frame with lon lat coordinates columns
#' @param cellsize positive numeric size of raster cells in meters
#' @param bandwidth positive numeric bandwidth in units of x
#' @param env_data SpatRaster object of environmental data
#' @param scale_01 boolean argument if activity data should be rescaled to 0-1 values range
#' @param data_extent TODO
#' @param input_crs character or terra crs object specifying coordinate reference system of coordinates in x data frame
#' @param output_crs character or terra crs object of coordinate reference system of output
#' @param stats vector of characters statistics to be calculated. See terra::global. "count", "range" and "area" additionally added.
#'
#' @return list of SpatRaster result and data frame of statistics
#'
#' @examples
#'
#' statistics = c("count", "area", "min", "max", "range", "mean", "std", 'sum')
#'
#'
#' # activity space
#' exposure_KDE(x = geolife_sandiego, cellsize = 50, bandwidth = 25,
#'  scale_01 = TRUE, input_crs = "EPSG:4326", output_crs = "EPSG:32611",
#'  stats = statistics)
#'
#' #environmental exposure
#'
#' ndvi_data = terra::rast(system.file("extdata/landsat_ndvi.tif", package = "twsagps"))
#'
#' exposure_KDE(x = geolife_sandiego, cellsize = 50, bandwidth = 25,
#'  env_data = ndvi_data, scale_01 = TRUE, input_crs = "EPSG:4326",
#'  output_crs = "EPSG:32611", stats = statistics)
#'
#'
#'
#' @export
exposure_KDE = function(x, cellsize=NULL, bandwidth = 200, env_data=NULL,
                    scale_01 = FALSE, data_extent = NULL, # TODO extent
                    input_crs = "EPSG:4326", output_crs=NULL, stats=NULL){

  buff_const = 0.05 #percent of extent as bbox

  x_proj = start_processing(x,  env_data, data_extent, input_crs, output_crs)

  if (!is.null(env_data)){ # change env_data crs beforehand
    env_data_proj = terra::project(env_data, terra::crs(x_proj))
  }



  # get extent
  extent = terra::ext(x_proj)

  # buffer around extent
  x_const_kde = abs(terra::xmax(extent) - terra::xmin(extent)) * buff_const
  y_const_kde = abs(terra::ymax(extent) - terra::ymin(extent)) * buff_const

  # new extent
  new_extent = c(terra::xmin(extent) - x_const_kde, terra::xmax(extent) + x_const_kde,
                 terra::ymin(extent) - y_const_kde, terra::ymax(extent) + y_const_kde)
  #TODO small cellsize disproportion


  # implement cellsize > 0 condition
  if (is.numeric(cellsize)) { # cellsize included

    if (terra::is.lonlat(x_proj)){ # crs units in degrees
      warning("Cellsize is not stable - cells are not rectangular")

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


  # coords of every cell
  grid_coords = terra::crds(grid_rast)

  # point coords
  coords = terra::crds(x_proj)

  # when in degrees smoothing parameter not sqrt until resolving the output of kde
  if (terra::linearUnits(x_proj) == 0){
    kde_data = TDA::kde(coords, Grid = grid_coords, h = bandwidth)
  } else {
    # kde with sqrt bandwidth
    kde_data = TDA::kde(coords, Grid = grid_coords, h = sqrt(bandwidth))
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

  if (scale_01 == TRUE){ # scale_01 kde values

    rast_max = terra::minmax(kde_rast)[2,] # max
    # calculate normalization to 0-1 range
    kde_rast = (kde_rast)/(rast_max)
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
