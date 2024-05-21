#' Line Segment Method exposure
#' @description Line Segment method activity space and environmental exposure. Based on code from Micheal Garber's repository michaeldgarber/microclim-static-v-dynam.
#' In order to receive activity space ignore env_data argument.
#'
#' @param data Data.frame, SpatVector points or sf data.frame containing only POINTS.
#' @param bandwidth Positive numeric. Bandwidth in meters if output has a longtitude/latitude CRS or in the units of output CRS.
#' @param cellsize Positive numeric. Size of raster cells in meters if output has a longtitude/latitude CRS or in the units of output CRS.
#' @param time_data Data masking. Column of data containing POSIXct class objects.
#' @param env_data Stars, SpatRaster, SpatVector or sf. Spatial environmental data. Activity space is calculated when not set. When argument is a SpatVector or sf object vector data is rasterized to output raster using "sum" function.
#' @param output_crs Character or terra crs object. Coordinate Reference System of output. If not set and env_data is a SpatRaster env_data's CRS is used.
#' @param grid_extent Stars, SpatRaster, SpatExtent, sf bbox object or numeric vector of 4 length c(xmin, xmax, ymin, ymax). If stars or SpatRaster grid_extent is output's grid and cellsize argument is ignored. If SpatExtent, sf bbox object or vector grid_extent is output's extent. If cellsize is set it is preserved at the cost of extent.
#' @param normalize Boolean. If TRUE activity space SpatRaster is normalized
#' @param group_split Data masking. Column of data based on which it is grouped and split. For n groups output will be n SpatRasters.
#' @param filepath Character. Output filename. See terra::writeRaster().
#' @param norm_group Boolean. When normalize is TRUE, norm_method is "range" and group_split is set. If FALSE each SpatRaster is normalized seperately. If TRUE SpatRasters are normalized to SpatRaster with highest max value.
#' @param x Data masking. x or longitude coordinates column name if data is a data.frame.
#' @param y Data masking. y or latitude coordinates column name if data is a data.frame.
#' @param input_crs Character or terra crs object. Coordinate Reference System of data's coordinates if data is a data.frame.
#' @param norm_method Character. Normalization method. Four methods - "center", "scale", "standardize" and "range". Default is "range". See BBmisc::normalize().
#' @param NA_val Numeric. Value in x and y marked as NA if data is a data.frame.
#' @param time_unit Character. Unit of time weights of time_data. Ignored if time_data not specified. See difftime().
#' @param env_field Data masking. Column of env_data that env_data will be rasterized on. Ignored if env_data not a SpatVector or sf class.
#' @param env_buff Positive numeric. Optional buffer around SpatVector/sf env_data in meters if output has a longtitude/latitude CRS or in the units of the CRS. See terra::buffer(). Ignored if env_data not a SpatVector or sf class.
#'
#'
#' @return SpatRaster
#'
#' @examples
#'
#' # activity space
#' exposure_LS(data = geolife_sandiego, x = lon, y = lat, time_data = dateTime, time_unit = "mins",
#'   cellsize = 50, bandwidth = 200, normalize = TRUE, input_crs = "EPSG:4326",
#'   output_crs = "EPSG:32611")
#'
#' # split by date
#' exposure_LS(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, group_split = date,
#'    bandwidth = 200, normalize = TRUE, norm_group = FALSE, input_crs = "EPSG:4326",
#'    output_crs = "EPSG:32611")
#'
#' # split by date and define extent with time_only_norm
#' extent = c(478000, 484000, 3618000, 3627000)
#'
#' exposure_LS(data = geolife_sandiego, x = lon, y = lat, time_data = dateTime, time_unit = "mins",
#'   cellsize = 50, group_split = date, bandwidth = 200, normalize = TRUE, norm_group = TRUE,
#'   grid_extent = extent, input_crs = "EPSG:4326",
#'   output_crs = "EPSG:32611")
#'
#' #environmental exposure
#'
#' ndvi_data = terra::rast(system.file("extdata/landsat_ndvi.tif", package = "twiGPS"))
#'
#' exposure_LS(data = geolife_sandiego, x = lon, y = lat, time_data = dateTime, time_unit = "mins",
#'   cellsize = 50, bandwidth = 200, env_data = ndvi_data, normalize = FALSE,
#'   input_crs = "EPSG:4326", output_crs = "EPSG:32611")
#'
#' # environmental exposure - use rast grid and split by date with time_only_norm
#' exposure_LS(data = geolife_sandiego, x = lon, y = lat, time_data = dateTime, time_unit = "mins",
#'   group_split = date, bandwidth = 200, env_data = ndvi_data, normalize = TRUE,
#'   norm_group = TRUE, grid_extent = ndvi_data, input_crs = "EPSG:4326",
#'   output_crs = "EPSG:32611")
#' @export


exposure_LS = function(data, x, y, NA_val, time_data, time_unit = "mins",
                       cellsize, group_split, bandwidth, env_data, env_field, env_buff,
                       normalize = FALSE, norm_method = "range", norm_group = FALSE,
                       grid_extent, input_crs, output_crs, filepath){

  # handle stars objects rasters as env_data
  if (!missing(env_data)) {
    if (inherits(env_data, "stars")){
      env_data = terra::rast(env_data)
    } else if (inherits(env_data, "sf")) {
      env_data = terra::vect(env_data)
    } else if (!inherits(env_data, c("SpatVector", "SpatRaster"))){
      stop("Invalid env_data - env_data neither stars, sf, SpatVector nor SpatRaster class")
    }
  }

  #handle bandwidth
  if (missing(bandwidth)){
    stop("Missing bandwidth argument. Provide valid bandwidth")
  } else if (length(bandwidth) != 1 || !is.numeric(bandwidth) || bandwidth <= 0){
    stop("Invalid bandwidth argument - bandwidth is neither positive nor single numeric value")
  }

  # handle grid_extent
  if (!missing(grid_extent)){
    if (inherits(grid_extent, "stars")){
      grid_extent = terra::rast(grid_extent)
    } else if (inherits(grid_extent, "bbox") || (is.vector(grid_extent) &&
                                                 inherits(grid_extent, "numeric") && length(grid_extent) == 4)){
      grid_extent = terra::ext(grid_extent)
    } else if (!inherits(grid_extent, c("SpatExtent", "SpatRaster"))){
      stop("Invalid grid_extent - grid_extent neither stars, SpatRaster, SpatExtent, bbox class nor numeric vector of 4 length")
    }
  }

  # handle normalize and norm_group
  if (is.na(as.logical(normalize))){
    stop("Invalid normalize argument - normalize argument cannot be interpreted as boolean")
  }

  if (is.na(as.logical(norm_group))){
    stop("Invalid norm_group argument - norm_group argument cannot be interpreted as boolean")
  }

  # handle norm_method
  if (!norm_method %in% c("center", "scale", "standardize", "range") && normalize) {
    warning('Invalid norm_method - applying default normalization method "range"')
    norm_method = "range"
  }




  # get spatial data with correct crs
  if (!missing(data)){
    if (inherits(data, "data.frame")){
      if (all(!c(missing(x), missing(y)))){
        x_enq = rlang::quo_name(rlang::enquo(x))
        y_enq = rlang::quo_name(rlang::enquo(y))

        if (all(c(x_enq, y_enq) %in% colnames(data))){

          data_proj = start_processing(data = data, x = x_enq,
                                       y = y_enq, NA_val = NA_val,
                                       env_data = env_data, grid_extent = grid_extent,
                                       input_crs = input_crs, output_crs = output_crs)
        } else {
          stop("Invalid x or y arguments - x or y are not a column in data")
        }
      } else {
        stop('Missing x or y arguments for data "data.frame" class')
      }
    } else {
      data_proj = start_processing(data = data, env_data = env_data,
                                   grid_extent = grid_extent, input_crs = input_crs,
                                   output_crs = output_crs)
    }
  } else {
    stop("Missing data argument - provide valid data argument")
  }

  # if time_data remove points with NA time_data
  if (!missing(time_data)){
    time_cname = rlang::quo_name(rlang::enquo(time_data))
    if (!time_cname %in% terra::names(data_proj)){
      stop("Invalid time_data argument - time_data is not a column name in data")
    } else if (any(is.na(data_proj[time_cname]))){
      n_row = nrow(data_proj)
      data_proj = tidyterra::drop_na(data_proj,  time_cname)
      message(paste0("Removing points with NA ", time_cname, " values - removing ", n_row - nrow(data_proj), " rows"))
    }
  }

  if (!missing(grid_extent) && inherits(grid_extent, "SpatRaster")){
    if (terra::crs(data_proj) != terra::crs(grid_extent)){
      grid_extent = terra::project(grid_extent, terra::crs(data_proj))
    }
    grid_rast = grid_extent
  } else {
    grid_rast = calc_grid(x = data_proj, cellsize = cellsize,
                          env_data = env_data, grid_extent = grid_extent)
  }

  # if env_data is vector data - create optional buffer and rasterize to grid raster
  if (!missing(env_data) && inherits(env_data, "SpatVector")) {

    env_data = env_vect(env = env_data, env_buff = env_buff,
                        env_field = env_field, grid = grid_rast)
  }


  if (!missing(env_data)){ # env_data included
    # project env_data to grid_rast cellsize
    env_data_proj = terra::project(env_data, grid_rast)

    # replace vals of grid with env values
    #env_resamp = terra::resample(env_data_proj, grid_rast)
    terra::values(grid_rast) = terra::values(env_data_proj)
  }


  if (missing(group_split)) {
    data_iter = list(data_proj) # only one item for for loop
  } else {
    enq_group_split = rlang::quo_name(rlang::enquo(group_split))
    if (enq_group_split %in% terra::names(data_proj)){
      data_iter = terra::split(data_proj, enq_group_split) # split data_proj by group_split
      message(paste0("Data split by group into ", length(data_iter), " items"))
    } else {
      data_iter = list(data_proj)
      warning("Invalid group_split argument - group_split is not a column in data. Data not split")
    }

  }


  buff_svc = terra::vect()
  output = list()



  for (data_i in data_iter) {

    # create line segments from spatial points
    trajectories = terra::vect(trajectories_fun(data_i))

    # create buffers over line segments
    traj_buff = trajectories |>
      tidyterra::select(line_id) |>
      terra::buffer(bandwidth)


    traj_buff$act = 1


    if (!missing(time_data)){
      duration_line_id = data_i |>
        dplyr::mutate(time_elapsed = as.numeric(difftime(dplyr::lead(.data[[time_cname]]),
                                                         .data[[time_cname]], units = time_unit)),
                      line_id = dplyr::row_number()) |>
        dplyr::select(line_id, time_elapsed)
      # last NA value set to 0 for normalization to not set lowest time value to 0
      duration_line_id$time_elapsed[nrow(duration_line_id)] = 0

      #1st type of normalize - normalize time_elapsed without spatial reference
      if(normalize && (!norm_group || norm_method != "range" || length(data_iter) == 1)) {

        if (norm_method != "range" && norm_group) {
          message(paste0('Norm_method is "', norm_method, '" - norm_group is TRUE is applicable only for norm_method "range". Norm group argument ignored. Normalizing each group seperately'))
        }


        duration_line_id$time_elapsed = normalization(duration_line_id$time_elapsed,
                                                          method = norm_method)

      }
      duration_df = as.data.frame(duration_line_id)

      traj_buff = dplyr::left_join(traj_buff, duration_df, by = 'line_id')
      traj_buff = traj_buff[,-2]
      names(traj_buff)[2] = "act"

    }

    buff_svc = c(buff_svc, traj_buff)

  }
  #remove empty vector
  buff_svc = buff_svc[-1]

  if(normalize && !missing(time_data) && norm_group && length(data_iter) > 1 &&
     norm_method == "range"){

    if (length(data_iter) > 1){

      max_buff = terra::vect(buff_svc)$act |> max(na.rm = TRUE)

      for (i in 1:length(buff_svc)){
        buff_svc[i]$act = buff_svc[i]$act / max_buff
      }

      # buff_list = unlist(vapply(buff_list, function(x) {
      #   if (!(all(is.na(x$act)))) {
      #     x$act = x$act / max_buff
      #   }
      #   return(list(x))
      # }, c(new("SpatVectorCollection"))))
    }
  }

  if (inherits(buff_svc, "SpatVector")){
    buff_svc = c(buff_svc)
  }

  for (buff_n in  1:length(buff_svc)) {

    if (!missing(env_data)){ #all of computation necessary only when calculating env_data
      # extract values and weights (area overlap) from grid for each cell of buffer
      exac_extr = exactextractr::exact_extract(grid_rast, sf::st_as_sf(buff_svc[buff_n]), progress = FALSE)

      traj_extract = Map(function(x, id) {x$line_id = id
      return(x)}, x = exac_extr, seq_along(exac_extr)) |> dplyr::bind_rows() |>
        dplyr::as_tibble() |>
        dplyr::relocate(line_id) |>
        dplyr::rename(e = 2)


      # calculate summarised exposure for each buffer
      traj_extract_line_id = traj_extract |>
        dplyr::group_by(line_id) |>
        dplyr::summarise(
          #Jan 9, 2024 use R's built-in weighted.mean() function
          #instead of calculating weighted average manually
          end_weights=stats::weighted.mean(
            x=e,
            w=coverage_fraction,
            na.rm=TRUE),
          #These weights are based on the areal overlap, not time
          sum_weights = sum(coverage_fraction,na.rm=TRUE),
          n_pixel = dplyr::n() # number of observations corresponds to number of pixels per line segment
        ) |>
        dplyr::ungroup()

      # join weights with spatial buffer
      weight_buff = dplyr::inner_join(buff_svc[buff_n], traj_extract_line_id, by = 'line_id')

      weight_buff$end_weights = weight_buff$end_weights * weight_buff$act

      # rasterize results
      rast_segment = terra::rasterize(weight_buff, grid_rast, field = "end_weights", fun = "sum")
    } else {
      rast_segment = terra::rasterize(buff_svc[buff_n], grid_rast, field = "act", fun = "sum")
    }
    output = suppressWarnings(append(output, (rast_segment)))

  }

  if (!missing(filepath)) { # save raster
    terra::writeRaster(output, filename = filepath)
    message(paste0("Saving output to ", filepath))
  }


  return(output)

}
