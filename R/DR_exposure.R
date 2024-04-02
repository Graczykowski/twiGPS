#' Density Ranking exposure
#'
#' @description
#' Density Ranking method activity space and environmental exposure. Using modifed DR function from Yen-Chi Chen's repository density_ranking.
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
#'
#' @return list of SpatRaster result and data frame of statistics
#'
#' @examples
#'
#' statistics = c("count", "area", "min", "max", "range", "mean", "std", 'sum')
#'
#'
#' # activity space
#' DR_exposure(x = geolife_sandiego, cellsize = 50, bandwidth = 70,
#'  start_crs = "WGS84", end_crs = "EPSG:32611", stats = statistics)
#'
#' #environmental exposure
#' ndvi_data = terra::rast(system.file("extdata/landsat_ndvi.tif", package = "twsagps"))
#'
#' DR_exposure(x = geolife_sandiego, cellsize = 50, bandwidth = 70,
#'  env_data = ndvi_data, start_crs = "WGS84",
#'  end_crs = "EPSG:32611", stats = statistics)
#'
#' @export


DR_exposure = function(x, cellsize = NULL, bandwidth = 200, env_data=NULL,
                       normalize = FALSE, data_extent = NULL, # TODO extent
                       start_crs = "WGS84", end_crs=NULL, stats=NULL,
                       act_and_env=FALSE){ # TODO act_and_env

  buff_const = 0.05 #percent of extent as bbox

  x_proj = start_processing(x, env_data, data_extent, start_crs, end_crs)

  if (!is.null(env_data)){ # change env_data crs beforehand
    env_data_proj = terra::project(env_data, terra::crs(x_proj))
  }


  # get extent
  extent = terra::ext(x_proj)

  # buffer around extent
  x_const_dr = (terra::xmax(extent) - terra::xmin(extent)) * buff_const
  y_const_dr = (terra::ymax(extent) - terra::ymin(extent)) * buff_const

  # new extent
  new_extent = c(terra::xmin(extent) - x_const_dr, terra::xmax(extent) + x_const_dr,
                 terra::ymin(extent) - y_const_dr, terra::ymax(extent) + y_const_dr)



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




  # point coords
  coords = terra::geom(x_proj)[,3:4]

  # number of cells

  x_cells = terra::ncol(grid_rast)
  y_cells = terra::nrow(grid_rast)



  # calculate DR data
  dr_data = DR_simple(coords, kernel= "Gaussian", h = bandwidth,
                      x_res=x_cells, y_res=y_cells, xlim=new_extent[1:2],
                      ylim=new_extent[3:4])

  ### PROBABLY DOESNT MATTER
  # # creating empty raster to insert KDE values
  # len_lon = length(dr_data$x_grid)
  # len_lat = length(dr_data$y_grid)
  #
  # lon_min = min(dr_data$x_grid)
  # lon_max = max(dr_data$x_grid)
  # lat_min = min(dr_data$y_grid)
  # lat_max = max(dr_data$y_grid)
  #
  # dr_rast = terra::rast(nrows=len_lat, ncols=len_lon, xmin=lon_min, xmax=lon_max,
  #                       ymin=lat_min, ymax=lat_max, crs = terra::crs(x_proj))
  ###


  dr_rast = grid_rast
  # insert dr values to raster
  terra::values(dr_rast) = dr_data$gr_alpha

  if (normalize == TRUE){ # normalize dr values

    rast_minmax = terra::minmax(dr_rast) # minmax
    # calculate normalization to 0-1 range
    dr_rast = (dr_rast - rast_minmax[1,])/(rast_minmax[2,]-rast_minmax[1,])
  }

  dr_rast = terra::subst(dr_rast, from = 0, to = NA)
  dr_rast = terra::flip(dr_rast, direction='vertical') # flip raster


  if (!is.null(env_data)){ # calculate exposure
    # only for one point data check
    # env_max = terra::minmax(env_data)[2]
    # env_data_proj[env_data_proj > env_max] = NA
    ##
    env_data_resamp = terra::resample(env_data_proj, dr_rast)
    ##
    # env_data_resamp[env_data_resamp > env_max] = NA
    ##
    rast_env_dr = dr_rast * env_data_resamp

    # calculate env output
    output = output_calc(rast_env_dr, stats = stats)
  } else {
    # calculate activity output
    output = output_calc(dr_rast, stats = stats)

  }

  # generate output





  return(output)




}
