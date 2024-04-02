#' Point Overlay exposure
#'
#' @description
#' Point Overlay method activity space and environmental exposure. In order to receive activity space ignore env_data argument.
#'
#'
#' @param x data frame with lon lat coordinates columns
#' @param cellsize positive numeric size of raster cell in meters
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
#'
#' # activity space
#' PO_exposure(x = geolife_sandiego, cellsize = 50, normalize = TRUE,
#'  start_crs = "WGS84", end_crs = "EPSG:32611", stats = statistics)
#'
#' #environmental exposure
#'
#' ndvi_data = terra::rast(system.file("extdata/landsat_ndvi.tif", package = "twsagps"))
#'
#' PO_exposure(x = geolife_sandiego, cellsize = 50, env_data = ndvi_data,
#'  normalize = TRUE, start_crs = "WGS84", end_crs = "EPSG:32611",
#'  stats = statistics)
#'
#'
#'
#'
#' @export


PO_exposure = function(x, cellsize=NULL, env_data=NULL, normalize = FALSE,
                       data_extent = NULL, # TODO extent
                       start_crs = "WGS84", end_crs=NULL, stats=NULL,
                       act_and_env=FALSE){ # TODO act_and_env


  x_proj = start_processing(x, env_data, data_extent, start_crs, end_crs)

  if (!is.null(env_data)){ # change env_data crs beforehand
    env_data_proj = terra::project(env_data, terra::crs(x_proj))
  }

  # rasterize

  x_proj$rast = 1




  # implement cellsize > 0 condition
  if (is.numeric(cellsize)) { # cellsize included

    if (terra::linearUnits(x_proj) == 0){ # crs units in degrees
      # PUT IT IN SEPERATE FUNCTION

      extent = terra::ext(x_proj)
      dist_lon = geosphere::distm(c(extent[1], extent[3]), c(extent[2], extent[3]),
                                  fun = geosphere::distHaversine)
      dist_lat = geosphere::distm(c(extent[1], extent[3]), c(extent[1], extent[4]),
                                  fun = geosphere::distHaversine)
      # number of cells
      x_cells = (dist_lon / cellsize) |> as.integer()
      y_cells = (dist_lat / cellsize) |> as.integer()

      ### PROBABLY DOESNT MATTER
      # x_seq = seq(extent[1], extent[2], length.out = x_cells)
      # y_seq = seq(extent[3], extent[4], length.out = y_cells)
      #
      #
      # # params for empty raster
      # len_x <- length(x_seq)
      # len_y <- length(y_seq)
      #
      # # x, y limits
      # x_min = min(x_seq)
      # x_max = max(x_seq)
      # y_min = min(y_seq)
      # y_max = max(y_seq)
      ###

      # empty_rast for units in degrees
      grid_rast = terra::rast(crs = terra::crs(x_proj), nrows=y_cells,
                              ncols=x_cells, extent = extent)

    } else {
      grid_rast = terra::rast(crs = terra::crs(x_proj), extent = terra::ext(x_proj),
                               resolution = cellsize)
    }



  } else if (!is.null(env_data)){ #if incorrect cellsize and env_data exists

    grid_rast = env_data_proj
    terra::ext(grid_rast) = terra::ext(x_proj)

  }

  rast_points = terra::rasterize(x_proj, grid_rast, field = "rast", fun = sum)

  if (normalize == TRUE){

    rast_points = terra::subst(rast_points, from = NA, to = 0) # proper range for normalization
    rast_minmax = terra::minmax(rast_points) # minmax
    # calculate normalization to 0-1 range
    rast_points = (rast_points - rast_minmax[1,])/(rast_minmax[2,]-rast_minmax[1,])
    rast_points = terra::subst(rast_points, from = 0, to = NA) # insert NA

  }



  if (!is.null(env_data)){ # calculate exposure
    # only for one point data check
    # env_max = terra::minmax(env_data)[2]
    # env_data_proj[env_data_proj > env_max] = NA
    ##
    env_data_resamp = terra::resample(env_data_proj, rast_points)
    ## only for one point data check
    # env_data_resamp[env_data_resamp > env_max] = NA
    ##
    rast_env_points = rast_points * env_data_resamp


    # calculate env output
    output = output_calc(rast_env_points, stats = stats)
  } else {
    # calculate activity output
    output = output_calc(rast_points, stats = stats)
  }


  return(output)
}
