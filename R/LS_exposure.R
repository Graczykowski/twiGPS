#' Line Segment Method exposure
#' @description Based on code from Micheal garber's repository microclim-static-v-dynam
#'
#' @param x data frame with lon lat coordinate columns
#' @param day string in date format compatible with date column in x
#' @param time_data name of column in x containing POSIXct data class
#' @param cellsize size of raster cell in meters
#' @param buff_dist size of buffer im meter
#' @param env_data SparRaster object of envirinmental data
#' @param data_extent TODO
#' @param start_crs coordinate system of coordinates in x data frame
#' @param end_crs coordinate system of output
#' @param stats statistics calculated
#' @param act_and_env TODO
#'
#' @return list of SpatRaster result and list of statistics
#'
#' @export


# TODO optimalise function
# TODO whend result in WGS84 fix cellsize to fit degrees
LS_exposure = function(x, day=NULL, time_data = NULL, cellsize=100, buff_dist = 200,
                       env_data=NULL, data_extent = NULL, # TODO extent
                       start_crs = "WGS84", end_crs=NULL, stats=NULL,
                       act_and_env=FALSE){ # TODO act_and_env

  # processing dplyr argument
  time_data_null = dplyr::enquo(time_data)

  # get spatial data with correct crs
  x_proj = start_processing(x, day, env_data, data_extent, start_crs, end_crs)

  # create line segments from spatial points
  trajectories = terra::vect(trajectories_fun(x_proj))

  # create buffers over line segments
  traj_buff = trajectories |>
    dplyr::select(line_id) |>
    terra::buffer(buff_dist) # takes a lot of time


  # create grid raster

  if (is.numeric(cellsize) & cellsize > 0) { # cellsize included
    grid_rast = terra::rast(crs = terra::crs(trajectories), extent = terra::ext(trajectories),
                            resolution = cellsize, vals = 1) # vals of grid with weight 1

    if (!is.null(env_data)){ # env_data included
      # replace vals of grid with env values
      env_proj = terra::project(env_data, grid_rast)
      env_resamp = terra::resample(env_proj, grid_rast)
      terra::values(grid_rast) = terra::values(env_resamp)

    }

  } else if (!is.null(env_data)){ #if incorrect cellsize and env_data exists
    # grid data as env grid
    grid_rast = terra::project(env_data, terra::crs(trajectories))

    terra::ext(grid_rast) = terra::ext(trajectories) # ext as line segments
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
      dplyr::mutate(time_elapsed = as.integer(dplyr::lead({{time_data}}) - {{time_data}}),
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
