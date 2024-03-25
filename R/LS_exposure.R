#' Line Segment Method exposure
#' @description Line Segment method activity space and environmental exposure. Based on code from Micheal garber's repository microclim-static-v-dynam.
#' In order to receive activity space ignore env_data argument.
#'
#' @param x data frame with lon lat coordinates columns
#' @param time_data name of column in x containing POSIXct data class
#' @param time_unit character indicating unit of weights of time_data. Ignored if time_data not specified. See difftime()
#' @param cellsize positive numeric size of raster cell in meters
#' @param bandwidth positive numeric size of segments in meters
#' @param env_data SpatRaster object of environmental data
#' @param data_extent TODO
#' @param start_crs character or terra crs object specifying coordinate reference system of coordinates in x data frame
#' @param end_crs character or terra crs object of coordinate reference system of output
#' @param stats vector of characters statistics to be calculated. See terra::global. "count", "range" and "area" additionally added.
#'
#'
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
#' LS_exposure(x = geolife_sandiego, time_data = dateTime, time_unit = "mins",
#'  cellsize = 50, bandwidth = 200, start_crs = "WGS84", end_crs = "EPSG:32611",
#'   stats = statistics)
#'
#' #environmental exposure
#' data("landsat_ndvi")
#' ndvi_data = terra::rast(landsat_ndvi)
#'
#' LS_exposure(x = geolife_sandiego, time_data = dateTime, time_unit = "mins",
#'  cellsize = 50, bandwidth = 200, env_data = ndvi_data, start_crs = "WGS84",
#'  end_crs = "EPSG:32611", stats = statistics)
#'
#'
#' @export


# TODO optimalise function
# TODO whend result in WGS84 fix cellsize to fit degrees
LS_exposure = function(x, time_data = NULL, time_unit = "mins",
                       cellsize=NULL, bandwidth = 200,
                       env_data=NULL, data_extent = NULL, # TODO extent
                       start_crs = "WGS84", end_crs=NULL, stats=NULL,
                       act_and_env=FALSE){ # TODO act_and_env

  # processing dplyr argument
    time_data_null = dplyr::enquo(time_data)

  # get spatial data with correct crs
  x_proj = start_processing(x, env_data, data_extent, start_crs, end_crs)

  if (!is.null(env_data)){ # change env_data crs beforehand
    env_data_proj = terra::project(env_data, terra::crs(x_proj))
  }

  # create line segments from spatial points
  trajectories = terra::vect(trajectories_fun(x_proj))

  # create buffers over line segments
  traj_buff = trajectories |>
    tidyterra::select(line_id) |>
    terra::buffer(bandwidth) # takes a lot of time

  #
  extent_buff = terra::ext(traj_buff)



  # create grid raster

  # implement cellsize > 0 condition
  if (is.numeric(cellsize)) { # cellsize included

    if (terra::linearUnits(x_proj) == 0){ # crs units in degrees


      dist_lon = geosphere::distm(c(extent_buff[1], extent_buff[3]), c(extent_buff[2], extent_buff[3]),
                                  fun = geosphere::distHaversine)
      dist_lat = geosphere::distm(c(extent_buff[1], extent_buff[3]), c(extent_buff[1], extent_buff[4]),
                                  fun = geosphere::distHaversine)
      # number of cells
      x_cells = (dist_lon / cellsize) |> as.integer()
      y_cells = (dist_lat / cellsize) |> as.integer()



      ### PROBABLY DOESNT MATTER
      # x_seq = seq(extent_buff[1], extent_buff[2], length.out = x_cells)
      # y_seq = seq(extent_buff[3], extent_buff[4], length.out = y_cells)


      # # params for empty raster
      # len_x <- length(x_seq)
      # len_y <- length(y_seq)

      # # x, y limits
      # x_min = min(x_seq)
      # x_max = max(x_seq)
      # y_min = min(y_seq)
      # y_max = max(y_seq)

      # #grid rast for units in degrees
      # grid_rast = terra::rast(crs = terra::crs(trajectories), nrows=y_cells,
      #                         ncols=x_cells, xmin=x_min, xmax=x_max, ymin=y_min,
      #                         ymax=y_max, vals = 1) # vals of grid with weight 1
      ### PROBABLY DOESNT MATTER

      #grid rast for units in degrees
      grid_rast = terra::rast(crs = terra::crs(trajectories), nrows=y_cells,
                              ncols=x_cells, extent = extent_buff,
                              vals = 1) # vals of grid with weight 1




    } else {
      grid_rast = terra::rast(crs = terra::crs(trajectories), extent = extent_buff,
                              resolution = cellsize, vals = 1) # vals of grid with weight 1
    }

    if (!is.null(env_data)){ # env_data included
      # replace vals of grid with env values
      env_resamp = terra::resample(env_data_proj, grid_rast)
      terra::values(grid_rast) = terra::values(env_resamp)

    }

  } else if (!is.null(env_data)){ #if incorrect cellsize and env_data exists
    # grid data as env grid
    grid_rast = env_data_proj

    terra::ext(grid_rast) = extent_buff # ext as line segments
  }

  # extract values and weights (area overlap) from grid for each cell of buffer
  traj_extract= grid_rast |>
    terra::extract(
      traj_buff,
      na.rm=TRUE,
      weights = TRUE
    ) |>
    dplyr::as_tibble() |>
    dplyr::rename(
      line_id = ID,#rename this to line id
      e=2#second column is the exposure.
    )

  # calculate summarised exposure for each buffer
  traj_extract_line_id = traj_extract |>
    dplyr::group_by(line_id) |>
    dplyr::summarise(
      #Jan 9, 2024 use R's built-in weighted.mean() function
      #instead of calculating weighted average manually
      e=stats::weighted.mean(
        x=e,
        w=weight,
        na.rm=T),
      #These weights are based on the areal overlap, not time
      sum_weights = sum(weight,na.rm=T),
      n_pixel = dplyr::n() # number of observations corresponds to number of pixels per line segment
    ) |>
    dplyr::ungroup()

  if (!(rlang::quo_is_null(time_data_null))){

    duration_line_id = x |>
      dplyr::mutate(time_elapsed = as.numeric(difftime(dplyr::lead({{time_data}}),
                                                       {{time_data}}, units = time_unit )),
                    line_id = dplyr::row_number()) |>
      dplyr::select(line_id, time_elapsed)

    traj_extract_line_id = traj_extract_line_id |>
      #now link in the time weight
      dplyr::left_join(duration_line_id,by=c("line_id")) |>
      dplyr::mutate(
        #calculate a weight that considers both area overlap and time
        end_weights = e * time_elapsed
      )
  } else {
      traj_extract_line_id = traj_extract_line_id |>
        dplyr::mutate(end_weights = e) # end result is exposure without time
    }
  # join weights with spatial buffer
  weight_buff = dplyr::inner_join(traj_buff, traj_extract_line_id, by = 'line_id')

  # rasterize results
  rast_segment =terra::rasterize(weight_buff, grid_rast, field = "end_weights", fun = sum)

  # generate output
  output = output_calc(rast_segment, stats = stats)

  return(output)

}
