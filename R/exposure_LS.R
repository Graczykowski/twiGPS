#' Line Segment Method exposure
#' @description Line Segment method activity space and environmental exposure. Based on code from Micheal Garber's repository michaeldgarber/microclim-static-v-dynam.
#' In order to receive activity space ignore env_data argument.
#'
#' @param data Data.frame, SpatVector points or sf data.frame containing only POINTS.
#' @param x Data masking. x or lontitude coordinates column if data is a data.frame.
#' @param y Data masking. y or latitude coordinates column if data is a data.frame.
#' @param NA_val Numeric. Value in x and y marked as NA if data is a data.frame.
#' @param time_data Data masking. Column of data containing POSIXct class objects.
#' @param time_unit Character. Unit of time weights of time_data. Ignored if time_data not specified. See difftime().
#' @param cellsize Positive numeric. Size of raster cells in meters if output has a longtitude/latitude CRS or in the units of output CRS.
#' @param group_split Data masking. Column of data based on which it is grouped and split. For n groups output will be n SpatRasters.
#' @param bandwidth Positive numeric. Bandwidth in units of output Coordinate Reference System.
#' @param env_data Stars, SpatRaster, SpatVector or sf. Spatial environmental data. Activity space is calculated when not set. When argument is a SpatVector or sf object vector data is rasterized to output raster using "sum" function.
#' @param env_field Data masking. Column of env_data that env_data will be rasterized on. Ignored if env_data not a SpatVector or sf class.
#' @param env_buff Positive numeric. Optional buffer around SpatVector/sf env_data in meters if output has a longtitude/latitude CRS or in the units of the CRS. See terra::buffer(). Ignored if env_data not a SpatVector or sf class.
#' @param normalize Boolean. If activity data should be normalized.
#' @param norm_method Character. Normalization method. Four methods - "center", "scale", "standardize" and "range". Default is "range". See BBmisc::normalize().
#' @param norm_group Boolean. When normalize is TRUE, norm_method is "range" and group_split is set. If FALSE each SpatRaster is rescaled seperately. If TRUE SpatRasters are rescaled to SpatRaster with highest max value.
#' @param time_only_norm Boolean. Time type of normalization. If TRUE only time_data is normalized, normalized raster values may be imprecise. If FALSE and norm_method is "range" normalization is conducted with regard to spatial distribution, normalized raster values are accurate. TEMPORARY
#' @param grid_extent Stars, SpatRaster, SpatExtent, sf bbox object or numeric vector of 4 length c(xmin, xmax, ymin, ymax). If stars or SpatRaster grid_extent is output's grid and cellsize argument is ignored. If SpatExtent, sf bbox object or vector grid_extent is output's extent. If cellsize is set it is preserved at the cost of extent.
#' @param input_crs Character or terra crs object. Coordinate Reference System of data's coordinates if data is a data.frame.
#' @param output_crs Character or terra crs object. Coordinate Reference System of output. If not set and env_data is a SpatRaster env_data's CRS is used.
#' @param filepath Character. Output filename. See terra::writeRaster().
#'
#'
#' @return SpatRaster
#'
#' @examples
#'
#' # activity space
#' exposure_LS(data = geolife_sandiego, x = lon, y = lat, time_data = dateTime, time_unit = "mins",
#'   cellsize = 50, bandwidth = 200, normalize = TRUE, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
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
#'   time_only_norm = TRUE, grid_extent = extent, input_crs = "EPSG:4326",
#'   output_crs = "EPSG:32611")
#'
#' #environmental exposure
#'
#' ndvi_data = terra::rast(system.file("extdata/landsat_ndvi.tif", package = "twsagps"))
#'
#' exposure_LS(data = geolife_sandiego, x = lon, y = lat, time_data = dateTime, time_unit = "mins",
#'   cellsize = 50, bandwidth = 200, env_data = ndvi_data, normalize = FALSE,
#'   input_crs = "EPSG:4326", output_crs = "EPSG:32611")
#'
#' # environmental exposure - use rast grid and split by date with time_only_norm
#' exposure_LS(data = geolife_sandiego, x = lon, y = lat, time_data = dateTime, time_unit = "mins",
#'   group_split = date, bandwidth = 200, env_data = ndvi_data, normalize = TRUE,
#'   norm_group = TRUE, time_only_norm = TRUE, grid_extent = ndvi_data, input_crs = "EPSG:4326",
#'   output_crs = "EPSG:32611")
#' @export


# TODO optimalise function
exposure_LS = function(data, x, y, NA_val, time_data, time_unit = "mins",
                       cellsize, group_split, bandwidth = 200, env_data, env_field, env_buff,
                       normalize = FALSE, norm_method = "range", norm_group = FALSE, time_only_norm = FALSE,
                       grid_extent, input_crs, output_crs, filepath){

  # handle stars objects rasters as env_data
  if (!missing(env_data)) {
    if (any(class(env_data) == "stars")){
      env_data = terra::rast(env_data)
    } else if (any(class(env_data) == "sf")) {
      env_data = terra::vect(env_data)
    } else if (any(!class(env_data) %in% c("SpatVector", "SpatRaster"))){
      stop("Invalid env_data - env_data neither stars, sf, SpatVector nor SpatRaster class")
    }
  }

  # handle grid_extent
  if (!missing(grid_extent)){
    if (any(class(grid_extent) == "stars")){
      grid_extent = terra::rast(grid_extent)
    } else if (any(class(grid_extent) == "bbox" || (is.vector(grid_extent) &&
                                                    all(class(grid_extent) == "numeric") && length(grid_extent) == 4))){
      grid_extent = terra::ext(grid_extent)
    } else if (any(!class(grid_extent) %in% c("SpatExtent", "SpatRaster"))){
      stop("Invalid grid_extent - grid_extent neither stars, SpatRaster, SpatExtent, bbox class nor numeric vector of 4 length")
    }
  }

  # handle norm_method
  if (!norm_method %in% c("center", "scale", "standardize", "range") && normalize == TRUE) {
    warning('Ivnalid norm_method. Applying default normalization method "range"')
    norm_method = "range"
  }

  # get spatial data with correct crs
  if (all(class(data) == "data.frame")){
    if (all(!c(missing(x), missing(y)))){
      x_enq = rlang::enquo(x)
      y_enq = rlang::enquo(y)

      data_proj = start_processing(data = data, x = rlang::quo_name(x_enq),
                                   y = rlang::quo_name(y_enq), NA_val = NA_val,
                                   env_data = env_data, grid_extent = grid_extent,
                                   input_crs = input_crs, output_crs = output_crs)
    } else {
      stop('Incorrect x and y arguments for data "data.frame" class')
    }
  } else {
    data_proj = start_processing(data = data, env_data = env_data,
                                 grid_extent = grid_extent, input_crs = input_crs,
                                 output_crs = output_crs)
  }

  # if time_data remove points with NA time_data
  if (!missing(time_data)){
    enq_time = rlang::enquo(time_data)
    time_cname = rlang::quo_name(enq_time)
    if (any(is.na(data_proj[time_cname]))){
      message(paste0("Removing points with NA ", time_cname, " values"))
      data_proj = tidyterra::drop_na(data_proj,  time_cname)
    }
  }

  if (!missing(grid_extent) && any(class(grid_extent) == "SpatRaster")){
    grid_rast = grid_extent
  } else {
    grid_rast = calc_grid(x = data_proj, cellsize = cellsize, env_data = env_data,
                          grid_extent = grid_extent, bandwidth = bandwidth, is_LS = TRUE)
  }


  # if env_data is vector data - create optional buffer and rasterize to grid raster
  if (!missing(env_data) && any(class(env_data) == "SpatVector")) {

    env_data = env_vect(env = env_data, env_buff = env_buff,
                        env_field = env_field, grid = grid_rast)
  }


  if (!missing(env_data) && any(class(env_data) == "SpatRaster")){ # env_data included
    # project env_data to grid_rast cellsize
    env_data_proj = terra::project(env_data, grid_rast)

    # replace vals of grid with env values
    #env_resamp = terra::resample(env_data_proj, grid_rast)
    terra::values(grid_rast) = terra::values(env_data_proj)
  }


  if (missing(group_split)) {
    data_iter = list(data_proj) # only one item for for loop
  } else {
    enq_group_split = rlang::enquo(group_split)
    data_iter = terra::split(data_proj, rlang::quo_name(enq_group_split)) # split data_proj by group_split
    message(paste0("Data split by group into ", length(data_iter), " items"))
  }


  if (length(data_iter) > 1 && !missing(filepath)) {
    # list output

    # for iterating file names
    file_ext = stringi::stri_extract(filepath, regex = "\\.(\\w+)$")
    file_no_ext = substr(filepath, 1, nchar(filepath) - nchar(file_ext))

    file_vect = rep(file_no_ext, length.out = length(data_iter))
    file_vect = paste0(file_vect, "_", seq_along(file_vect), file_ext)
  }

  output = list()
  buff_list = list()


  for (data_i in data_iter) {

    # create line segments from spatial points
    trajectories = terra::vect(trajectories_fun(data_i))

    # create buffers over line segments
    traj_buff = trajectories |>
      tidyterra::select(line_id) |>
      terra::buffer(bandwidth)


    traj_buff$act = 1


    if (!missing(time_data)){
      #time_data_null = dplyr::enquo(time_data)
      duration_line_id = data_i |>
        dplyr::mutate(time_elapsed = as.numeric(difftime(dplyr::lead(.data[[time_cname]]),
                                                         .data[[time_cname]], units = time_unit)),
                      line_id = dplyr::row_number()) |>
        dplyr::select(line_id, time_elapsed)


      #1st type of normalize - normalize time_elapsed without spatial reference
      if(normalize == TRUE && time_only_norm == TRUE &&
         (norm_group != TRUE || norm_method != "range" || missing(group_split))) {

        if (norm_method == "range") {
          message('Norm_method is "range". Norm_group argument ignored. Normalizing each group seperately')
        }

        duration_line_id$time_elapsed = BBmisc::normalize(duration_line_id$time_elapsed,
                                                          method = norm_method, margin = 1L)

      }
      duration_df = as.data.frame(duration_line_id)

      traj_buff = dplyr::left_join(traj_buff, duration_df, by = 'line_id')
      traj_buff = traj_buff[,-2]
      names(traj_buff)[2] = "act"

    }

    if(normalize == TRUE && (time_only_norm == FALSE || missing(time_data)) &&
       (norm_group == FALSE || missing(group_split))){
      #2nd type of normalize - normalize time_elapsed or buffers using spatial reference - rasterizing and reading min max data from raster

      if (norm_method != "range") {
        warning(paste0('Norm_method "', norm_method, '" is inapplicable to not only time normalization. Using norm_method "range"'))
      }

      # slightly varying output values (max very close to 1)

      norm_rast = terra::rasterize(traj_buff, grid_rast, field = "act", fun = "sum")
      #read max val
      max_norm = terra::minmax(norm_rast)[2]

      # normalize
      if (length(unique(traj_buff$act)) == 1){
        #if vector is constant adjust BBmisc
        traj_buff$act = BBmisc::normalize(traj_buff$act, method = "range",
                                          margin = 1L, range = c(0, max(traj_buff$act) / max_norm * 2))
      } else{
        traj_buff$act = BBmisc::normalize(traj_buff$act, method = "range",
                                          margin = 1L, range = c(0, max(traj_buff$act) / max_norm))
      }

    }

    buff_list[length(buff_list) + 1] = list(traj_buff)

  }

  if(normalize == TRUE && time_only_norm == TRUE && !missing(time_data) &&
     norm_group == TRUE && !missing(group_split) && norm_method == "range"){
    # 1st normalize variant for all rasts


    if (length(data_iter) > 1){

      max_buff = (do.call(rbind, buff_list))$act |> max(na.rm = TRUE)
      buff_list = sapply(buff_list, function(x) {
        if (!(all(is.na(x$act)))) {
          x$act = BBmisc::normalize(x$act, method = "range",
                                               margin = 2L,
                                               range = c(0, max(x$act, na.rm = TRUE) / max_buff))
        }
        return(x)
      })

    } else {

      max_buff = buff_list[[1]]$act |> max(na.rm = TRUE)
      buff_list = sapply(buff_list, function(x) {
        if (!(all(is.na(x$act)))) {
          x$act = BBmisc::normalize(x$act, method = "range",
                                               margin = 2L,
                                               range = c(0, max(x$act, na.rm = TRUE) / max_buff))
        }
        return(x)
      })
    }
  } else if (normalize == TRUE && (time_only_norm == FALSE || missing(time_data)) &&
             norm_group == TRUE && !missing(group_split)){
    # 2nd normalize variant for all rasts
    if (norm_method != "range") {
      warning(paste0('Norm_method "', norm_method, '" is inapplicable to not only time normalization. Using norm_method "range"'))
    }
    norm_list = sapply(buff_list, terra::rasterize, grid_rast, field = "act", fun = "sum")
    max_multiple_norm = max(sapply(norm_list, terra::minmax)[2,], na.rm = TRUE)
    buff_list = sapply(buff_list, function(x) {
      if (!(all(is.na(x$act)))) {
        x$act = BBmisc::normalize(x$act, method = "range",
                                  margin = 2L,
                                  range = c(0, max(x$act, na.rm = TRUE) / max_multiple_norm))
      }
      return(x)
    })
  }

  for (buff in buff_list) {

    if (!missing(env_data)){ #all of computation necessary only when calculating env_data
      # extract values and weights (area overlap) from grid for each cell of buffer
      traj_extract= grid_rast |>
        terra::extract(# takes a lot of time
          buff,
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
          end_weights=stats::weighted.mean(
            x=e,
            w=weight,
            na.rm=TRUE),
          #These weights are based on the areal overlap, not time
          sum_weights = sum(weight,na.rm=TRUE),
          n_pixel = dplyr::n() # number of observations corresponds to number of pixels per line segment
        ) |>
        dplyr::ungroup()

      # join weights with spatial buffer
      weight_buff = dplyr::inner_join(buff, traj_extract_line_id, by = 'line_id')

      weight_buff$end_weights = weight_buff$end_weights * weight_buff$act

      # rasterize results
      rast_segment = terra::rasterize(weight_buff, grid_rast, field = "end_weights", fun = "sum")
    } else {
      rast_segment = terra::rasterize(buff, grid_rast, field = "act", fun = "sum")
    }
    if (length(data_iter) == 1) { # if one group or no group output is not a list
      output = rast_segment
    } else {
      output[length(output) + 1] = as.list(rast_segment)
    }
  }

  if (!missing(filepath)) { # save raster
    if (length(data_iter) > 1){ # update filepath

      mapply(function(x, y) {
        terra::writeRaster(x, filename = y)
        message(paste0("Saving output to ", y))
      }, x = output, y = file_vect)

    } else {

      terra::writeRaster(end_rast, filename = filepath)
      message(paste0("Saving output to ", filepath))
    }
  }


  return(output)

}
