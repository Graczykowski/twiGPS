#' Density Ranking in Projected Coordinate Systems
#' @description
#' Code working fast with results to check, working well in projected coordinate systems, WGS84 yet to do
#' Using modifed DR function from Yen-Chi Chen's repository density_ranking
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





DR_proj = function(x, day=NULL, cellsize=100, bandwidth = 200, env_data=NULL,
                    normalize = FALSE, data_extent = NULL, # TODO extent
                    start_crs = "WGS84", end_crs=NULL, stats=NULL,
                    act_and_env=FALSE){ # TODO act_and_env

  buff_const = 0.05 #percent of extent as bbox

  x_proj = start_processing(x, day, env_data, data_extent, start_crs, end_crs)


  # get extent
  extent = terra::ext(x_proj)

  # buffer around extent
  x_const_dr = (terra::xmax(extent) - terra::xmin(extent)) * buff_const
  y_const_dr = (terra::ymax(extent) - terra::ymin(extent)) * buff_const

  # new extent
  new_extent = c(terra::xmin(extent) - x_const_dr, terra::xmax(extent) + x_const_dr,
                 terra::ymin(extent) - y_const_dr, terra::ymax(extent) + y_const_dr)

  if (!is.null(env_data)){ # change env_data crs beforehand
    env_data_proj = terra::project(env_data, terra::crs(x_proj))
  }



  if (is.numeric(cellsize) & cellsize > 0) { # cellsize included

    grid_rast = terra::rast(crs = terra::crs(x_proj))
    terra::ext(grid_rast) = new_extent # ext before cellsize to avoid cellsize disproportion
    terra::res(grid_rast) = cellsize
  } else if (!is.null(env_data)){ #if incorrect cellsize and env_data exists

    grid_rast = env_data_proj
    terra::ext(grid_rast) = new_extent # ext before cellsize to avoid cellsize disproportion
    terra::res(grid_rast) = terra::res(env_data_proj)
  }



  # coordinate seq
  x_seq = seq(new_extent[1], new_extent[2], length.out = terra::ncol(grid_rast))
  y_seq = seq(new_extent[3], new_extent[4], length.out = terra::nrow(grid_rast))

  # coords of every cell
  expand_grid = expand.grid(x_seq,y_seq)

  # point coords
  coords = terra::geom(x_proj)[,3:4]

  # number of cells
  x_cells = terra::ncol(grid_rast)
  y_cells = terra::nrow(grid_rast)

  # dr with sqrt bandwidth
  dr_data = DR_simple(coords, kernel= "Gaussian", h = bandwidth,
                      x_res=x_cells, y_res=y_cells, xlim=new_extent[1:2],
                      ylim=new_extent[3:4])




  # params for empty raster
  len_x <- length(dr_data$x_grid)
  len_y <- length(dr_data$y_grid)

  # x, y limits
  x_min = min(dr_data$x_grid)
  x_max = max(dr_data$x_grid)
  y_min = min(dr_data$y_grid)
  y_max = max(dr_data$y_grid)

  # empty rast for dr
  dr_rast = terra::rast(nrows=len_y, ncols=len_x, xmin=x_min, xmax=x_max, ymin=y_min,
                         ymax=y_max, crs = terra::crs(x_proj))


  if (normalize == TRUE){ # normalize dr values
    dr_vals = BBmisc::normalize(as.vector(dr_data$gr_alpha), method = "range", range = c(0, 1),
                                 margin = 1L, on.constant = "quiet")
  } else {
    dr_vals = dr_data$gr_alpha
  }

  # insert dr values to raster
  terra::values(dr_rast) = dr_vals
  dr_rast = terra::subst(dr_rast, from = 0, to = NA)
  dr_rast = terra::flip(dr_rast, direction='vertical') # flip raster



  if (!is.null(env_data)){ # calculate exposure
    env_data_resamp = terra::resample(env_data_proj, dr_rast)
    dr_env_rast = dr_rast * env_data_resamp


    # calculate env output
    output = output_calc(dr_rast, env_rast = dr_env_rast, stats = stats)
  } else {
    # calculate activity output
    output = output_calc(dr_rast, stats = stats)
  }

  return(output)

}
