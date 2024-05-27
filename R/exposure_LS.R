#' Calculate Line Segment exposure of point spatial data
#'
#' @description
#' `exposure_LS()` calculates activity space and environmental exposure using Line Segment method.
#' Line Segment function is based on code from Micheal Garber's repository michaeldgarber/microclim-static-v-dynam.
#'
#' This method creates trajectories following points of spatial data and creates buffers around each trajectory.
#' Subsequently, values for every buffer are extracted from `env_data`.
#' Finally, `time_data` weights are incorporated and converted to SpatRaster with [terra::rasterize()] by `fun = "sum"`.
#' Setting `env_data` argument determines whether environmental exposure or activity space is computed.
#'
#'
#' For calculating statistics of output raster use [exposure_stats()]
#'
#'
#' @param data Data.frame, SpatVector points or sf data.frame containing only POINTS.
#' @param coords Character. Vector of column names of x and y coordinates if `data` is a data.frame.
#' @param bandwidth Positive numeric. Bandwidth in meters if output has a longtitude/latitude CRS or in the units of output's CRS.
#' @param cellsize Positive numeric. Size of raster cells in meters if output has a longtitude/latitude CRS or in the units of output CRS.
#' @param time_data Data masking. Column of data containing POSIXct class objects.
#' @param env_data Stars, SpatRaster, SpatVector or sf. Spatial environmental data. Activity space is calculated when not set. When argument is a SpatVector or sf object vector data is rasterized to output raster using `"sum"` function.
#' @param output_crs Character or terra crs object. CRS of output. If not set and `env_data` is a SpatRaster `env_data`'s CRS is used.
#' @param input_crs Character or terra crs object. CRS of `data`'s coordinates if `data` is a data.frame.
#' @param grid_extent Stars, SpatRaster, SpatExtent, sf bbox object or numeric vector of 4 length `c(xmin, xmax, ymin, ymax)`. If stars or SpatRaster `grid_extent` is output's grid and `cellsize` argument is ignored. If SpatExtent, sf bbox object or vector `grid_extent` is output's extent. Then, if `cellsize` is set it is preserved at the cost of extent.
#' @param normalize Character. If set activity space SpatRaster is normalized with specified method. Four methods - "center", "scale", "standardize" and "range". See [BBmisc::normalize()].
#' @param group_split Data masking. Column of data based on which it is grouped and split. For n groups output will be SpatRasters with n layers.
#' @param norm_group Boolean. Applicable only When normalize is "range" and group_split is set. If FALSE each layer of SpatRaster is normalized seperately. If TRUE all layer are normalized to layer with highest max value (default FALSE).
#' @param filepath Character. Output filename. See [terra::writeRaster()].
#' @param NA_val Numeric. Value in x and y marked as NA if `data` is a data.frame.
#' @param time_unit Character. Unit of time weights of `time_data` (default `"mins"`). Ignored if `time_data` not specified. See [difftime()].
#' @param env_field Data masking. Column of `env_data` that `env_data` will be rasterized on. Ignored if `env_data` not a SpatVector or sf class.
#' @param env_buff Positive numeric. Optional buffer around SpatVector/sf `env_data` in meters if output has a longtitude/latitude CRS or in the units of the CRS. Ignored if `env_data` not a SpatVector or sf class. See [terra::buffer()].
#' @param verbose Boolean. If FALSE amount of output is reduced (default TRUE).
#'
#' @return SpatRaster
#'
#' @references Jankowska, Marta M., Jiue-An Yang, Nana Luo, Chad Spoon, and Tarik Benmarhnia. 2023. “Accounting for Space, Time, and Behavior Using GPS Derived Dynamic Measures of Environmental Exposure.” Health & Place 79 (January): 102706. https://doi.org/10.1016/j.healthplace.2021.102706.
#'
#' @examples
#' # activity space for data.frame 'data' wihout time data
#' exposure_LS(data = geolife_sandiego, coords = c("lon", "lat"),
#'   bandwidth = 200, cellsize = 50, output_crs = "EPSG:32611",
#'   input_crs = "EPSG:4326")
#'
#' #SpatVector data
#' geolife_vect = terra::vect(geolife_sandiego,
#'                            geom = c("lon", "lat"),
#'                            crs = "EPSG:4326") |>
#'                            terra::project("EPSG:32611")
#'
#' # time_data
#' exposure_LS(data = geolife_vect, bandwidth = 200,
#'   cellsize = 50, time_data = dateTime)
#'
#' # normalize "range"
#' exposure_LS(data = geolife_vect, bandwidth = 200,
#'   cellsize = 50, time_data = dateTime, normalize = "range")
#'
#' # normalize "center"
#' exposure_LS(data = geolife_vect, bandwidth = 200,
#'   cellsize = 50, time_data = dateTime, normalize = "center")
#'
#' # split by date
#' exposure_LS(data = geolife_vect, bandwidth = 200,
#'   cellsize = 50, time_data = dateTime,
#'   normalize = "range", group_split = date)
#'
#' # split by date, normalize by group and define extent
#' extent = c(478000, 484000, 3618000, 3627000)
#'
#' exposure_LS(data = geolife_vect, bandwidth = 200,
#'   cellsize = 50, time_data = dateTime,
#'   grid_extent = extent, normalize = "range",
#'   group_split = date, norm_group = TRUE)
#'
#' # environmental exposure
#' ndvi_data = terra::rast(system.file("extdata/landsat_ndvi.tif",
#'                                      package = "twiGPS"))
#'
#' exposure_LS(data = geolife_vect, bandwidth = 200,
#'  cellsize = 50, time_data = dateTime, env_data = ndvi_data)
#'
#' # environmental exposure - use rast grid and split by date
#' exposure_LS(data = geolife_vect, bandwidth = 200,
#'   time_data = dateTime, env_data = ndvi_data,
#'   grid_extent = ndvi_data, normalize = "range",
#'   group_split = date)
#' @seealso [exposure_PO()], [exposure_KDE()], [exposure_DR()]
#' @export


exposure_LS = function(data, coords, bandwidth, cellsize, time_data, env_data,
                       output_crs, input_crs, grid_extent, normalize, group_split,
                       norm_group = FALSE, filepath, NA_val, time_unit = "mins",
                       env_field, env_buff, verbose = TRUE){


  # handle verbose
  if (is.na(as.logical(verbose))){
    stop("Argument 'verbose' should be interpretable as boolean")
  }

  # handle stars objects rasters as env_data
  if (!missing(env_data)) {
    if (inherits(env_data, "stars")){
      env_data = terra::rast(env_data)
    } else if (inherits(env_data, "sf")) {
      env_data = terra::vect(env_data)
    } else if (!inherits(env_data, c("SpatVector", "SpatRaster"))){
      stop("Argument 'env_data' should be stars, sf, SpatVector or SpatRaster class")
    }
  }

  if (missing(bandwidth)){
    stop("Argument 'bandwidth' is missing - valid 'bandwidth' argument should be provided")
  } else if (length(bandwidth) != 1 || !is.numeric(bandwidth) || bandwidth <= 0){
    stop("Argument 'bandwidth' should be positive and single numeric value")
  }

  # handle grid_extent
  if (!missing(grid_extent)){
    if (inherits(grid_extent, "stars")){
      grid_extent = terra::rast(grid_extent)
    } else if (inherits(grid_extent, "bbox") || (is.vector(grid_extent) &&
                                                 inherits(grid_extent, "numeric") && length(grid_extent) == 4)){
      grid_extent = terra::ext(grid_extent)
    } else if (!inherits(grid_extent, c("SpatExtent", "SpatRaster"))){
      stop("Argument 'grid_extent' should be stars, SpatRaster, SpatExtent, sf bbox class or numeric vector of 4 length")
    }
  }

  # handle normalize
  if (!missing(normalize) && !normalize %in% c("center", "scale", "standardize", "range")) {
    if (verbose) {
      warning('Argumnet \'normalize\' should be "center", "scale", "standardize" or "range" - applying default normalization method "range"')
    }
    normalize = "range"
  }

  # handle norm_group
  if (is.na(as.logical(norm_group))){
    stop("Argument 'norm_group' should be interpretable as boolean")
  }


  # get spatial data with correct crs
  if (!missing(data)){
    if (inherits(data, "data.frame")){
      if (!missing(coords)){

        if (all(coords %in% colnames(data))){

          data_proj = start_processing(data = data, coords = coords, NA_val = NA_val,
                                       env_data = env_data, grid_extent = grid_extent,
                                       input_crs = input_crs, output_crs = output_crs,
                                       verbose = verbose)
        } else {
          stop("Column name or names in 'coords' argument are not in 'data'")
        }
      } else {
        stop('Argument \'coords\' is missing for \'data\' data.frame class - valid \'coords\' argument should be provided')
      }
    } else {
      data_proj = start_processing(data = data, env_data = env_data,
                                   grid_extent = grid_extent, input_crs = input_crs,
                                   output_crs = output_crs, verbose = verbose)
    }
  } else {
    stop("Argument 'data' is missing - valid 'data' argument should be provided")
  }

  # if time_data remove points with NA time_data
  if (!missing(time_data)){
    time_cname = rlang::quo_name(rlang::enquo(time_data))
    if (!time_cname %in% terra::names(data_proj)){
      stop("Column name in 'time_data' argument is not in 'data'")
    } else if (any(is.na(data_proj[time_cname][[1]]))){
      n_row = nrow(data_proj)
      data_proj = terra::na.omit(data_proj, time_cname)
      if (verbose) {
      number_rows = n_row - nrow(data_proj)
      warning(paste0("Removed ", number_rows, ifelse(number_rows == 1, " row", " rows"),
                     " containing NA in 'time_data' column name"))
      }
    }
  }

  if (!missing(grid_extent) && inherits(grid_extent, "SpatRaster")){
    if (terra::crs(data_proj) != terra::crs(grid_extent)){
      grid_extent = terra::project(grid_extent, terra::crs(data_proj))
    }
    grid_rast = grid_extent
  } else {
    grid_rast = calc_grid(x = data_proj, cellsize = cellsize, env_data = env_data,
                          grid_extent = grid_extent, verbose = verbose)
  }

  # if env_data is vector data - create optional buffer and rasterize to grid raster
  if (!missing(env_data) && inherits(env_data, "SpatVector")) {

    env_data = env_vect(env = env_data, env_buff = env_buff, env_field = env_field,
                        grid = grid_rast, verbose = verbose)
  }


  if (!missing(env_data)){ # env_data included
    # project env_data to grid_rast cellsize
    env_data_proj = terra::project(env_data, grid_rast)

    # replace vals of grid with env values
    terra::values(grid_rast) = terra::values(env_data_proj)
  }


  if (missing(group_split)) {
    data_iter = list(data_proj) # only one item for for loop
  } else {
    enq_group_split = rlang::quo_name(rlang::enquo(group_split))
    if (enq_group_split %in% terra::names(data_proj)){
      data_iter = terra::split(data_proj, enq_group_split) # split data_proj by group_split
      if (verbose) {
        message(paste0("Data is split into ", length(data_iter), " groups"))
      }
    } else {
      data_iter = list(data_proj)
      if (verbose) {
        warning("Column name in 'group_split' argument is not in 'data'. 'Data' is not split")
      }
    }
  }


  buff_svc = terra::vect()
  output = list()



  for (data_i in data_iter) {

    # create line segments from spatial points
    trajectories = terra::vect(trajectories_fun(data_i))

    # create buffers over line segments
    traj_buff = trajectories["line_id"] |>
      terra::buffer(bandwidth)


    traj_buff$act = 1


    if (!missing(time_data)){
      duration_line_id = data_i |> as.data.frame() |>
        dplyr::mutate(time_elapsed = as.numeric(difftime(dplyr::lead(.data[[time_cname]]),
                                                         .data[[time_cname]], units = time_unit)),
                      line_id = dplyr::row_number()) |>
        dplyr::select(line_id, time_elapsed)
      # last NA value set to 0 for normalization to not set lowest time value to 0
      duration_line_id$time_elapsed[nrow(duration_line_id)] = 0

      if(!missing(normalize) && (!norm_group || normalize != "range" || length(data_iter) == 1)) {

        if (normalize != "range" && norm_group && verbose) {
          message(paste0('Normalization method is "', normalize, '" - \'norm_group\' = TRUE is applicable only for normalization method "range". \'Norm_group\' argument ignored, each group is normalized seperately'))
        }


        duration_line_id$time_elapsed = normalization(duration_line_id$time_elapsed,
                                                          method = normalize)

      }
      #might use tidyterra for optimisation

      line_id_df = as.data.frame(traj_buff)
      buff_df = dplyr::left_join(line_id_df, duration_line_id, by = 'line_id')

      traj_buff$act = buff_df$time_elapsed


    }

    buff_svc = c(buff_svc, traj_buff)

  }
  #remove empty vector
  buff_svc = buff_svc[-1]

  if(!missing(normalize) && !missing(time_data) && norm_group && length(data_iter) > 1 &&
     normalize == "range"){

    if (length(data_iter) > 1){

      max_buff = terra::vect(buff_svc)$act |> max(na.rm = TRUE)

      for (i in 1:length(buff_svc)){
        buff_svc[i]$act = buff_svc[i]$act / max_buff
      }
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
        dplyr::relocate(3) |>
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

      weight_buff = terra::merge(buff_svc[buff_n], traj_extract_line_id)

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
    if (verbose) {
      message(paste0("Saving output to ", filepath))
    }
  }


  return(output)

}
