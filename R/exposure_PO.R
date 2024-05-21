#' Point Overlay exposure
#'
#' @description
#' Point Overlay method activity space and environmental exposure. In order to receive activity space ignore env_data argument.
#'
#'
#' @param data Data.frame, SpatVector points or sf data.frame containing only POINTS.
#' @param x Data masking. x or longitude coordinates column name if data is a data.frame.
#' @param y Data masking. y or latitude coordinates column name if data is a data.frame.
#' @param NA_val Numeric. Value in x and y marked as NA if data is a data.frame.
#' @param cellsize Positive numeric. Size of raster cells in meters if output has a longtitude/latitude CRS or in the units of output CRS.
#' @param group_split Data masking. Column of data based on which it is grouped and split. For n groups output will be n SpatRasters.
#' @param env_data Stars, SpatRaster, SpatVector or sf. Spatial environmental data. Activity space is calculated when not set. When argument is a SpatVector or sf object vector data is rasterized to output raster using "sum" function.
#' @param env_field Data masking. Column of env_data that env_data will be rasterized on. Ignored if env_data not a SpatVector or sf class.
#' @param env_buff Positive numeric. Optional buffer around SpatVector/sf env_data in meters if output has a longtitude/latitude CRS or in the units of the CRS. See terra::buffer(). Ignored if env_data not a SpatVector or sf class.
#' @param normalize Boolean. If TRUE activity space SpatRaster is normalized
#' @param norm_method Character. Normalization method. Four methods - "center", "scale", "standardize" and "range". Default is "range". See BBmisc::normalize().
#' @param norm_group Boolean. When normalize is TRUE, norm_method is "range" and group_split is set. If FALSE each SpatRaster is normalized seperately. If TRUE SpatRasters are normalized to SpatRaster with highest max value.
#' @param grid_extent Stars, SpatRaster, SpatExtent, sf bbox object or numeric vector of 4 length c(xmin, xmax, ymin, ymax). If stars or SpatRaster grid_extent is output's grid and cellsize argument is ignored. If SpatExtent, sf bbox object or vector grid_extent is output's extent. If cellsize is set it is preserved at the cost of extent.
#' @param input_crs Character or terra crs object. Coordinate Reference System of data's coordinates if data is a data.frame.
#' @param output_crs Character or terra crs object. Coordinate Reference System of output. If not set and env_data is a SpatRaster env_data's CRS is used.
#' @param filepath Character. Output filename. See terra::writeRaster().
#'
#' @return SpatRaster
#'
#' @examples
#'
#' # activity space
#' exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE,
#'  input_crs = "EPSG:4326", output_crs = "EPSG:32611")
#'
#' # split by date
#' exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, group_split = date,
#'    normalize = TRUE, norm_group = FALSE, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
#'
#' # split by date and define extent
#' extent = c(478000, 484000, 3618000, 3627000)
#'
#' exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, group_split = date,
#'   normalize = TRUE,norm_group = TRUE, grid_extent = extent, input_crs = "EPSG:4326",
#'   output_crs = "EPSG:32611")
#'
#' #environmental exposure
#'
#' ndvi_data = terra::rast(system.file("extdata/landsat_ndvi.tif", package = "twiGPS"))
#'
#' exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, env_data = ndvi_data,
#'  normalize = FALSE, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
#'
#' # environmental exposure - use rast grid and split by date
#' exposure_PO(data = geolife_sandiego, x = lon, y = lat, group_split = date,
#'   env_data = ndvi_data, normalize = FALSE, norm_group = TRUE, grid_extent = ndvi_data,
#'   input_crs = "EPSG:4326", output_crs = "EPSG:32611")
#'
#'
#'
#'
#' @export


exposure_PO = function(data, x, y, NA_val, cellsize, group_split, env_data, env_field, env_buff,
                       normalize = FALSE, norm_method = "range", norm_group = FALSE, grid_extent,
                       input_crs, output_crs, filepath){

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

  act_out = list()


  for (data_i in data_iter){
    # if each group should have seperate extent then output is a list rasts
    # if all groups should have same extent then output is rast with n layers

    ### UNCOMMENT IF EVERY RAST SHOULD HAVE SEPERATE EXTENT
    # if (missing(grid_extent)) {
    #   # crop ext of each rast
    #   grid_crop = terra::crop(grid_rast, terra::ext(data_i))
    #   rast_points = terra::rasterize(data_i, grid_crop,  fun = "length")
    # } else {
    #   rast_points = terra::rasterize(data_i, grid_rast,  fun = "length")
    # }


    ### UNCOMMENT IF EVERY RAST SHOULD HAVE SEPERATE EXTENT


    ### UNCOMMENT IF EVERY RAST SHOULD HAVE SAME EXTENT

    rast_points = terra::rasterize(data_i, grid_rast,  fun = "length")


    ### UNCOMMENT IF EVERY RAST SHOULD HAVE SAME EXTENT


    if (normalize && (!norm_group || norm_method != "range" || length(data_iter) == 1)){
      if (norm_method != "range" && norm_group) {
        message(paste0('Norm_method is "', norm_method, '" - norm_group is TRUE is applicable only for norm_method "range". Norm group argument ignored. Normalizing each group seperately'))
      }
      # calculate normalization to 0-1 range
      rast_points = normalization(rast_points, method = norm_method)

    }

    ### UNCOMMENT IF EVERY RAST SHOULD BE A SEPERATE ELEMENT IN LIST

    # act_out[length(act_out) + 1] = as.list(rast_points)

    ### UNCOMMENT IF EVERY RAST SHOULD BE A SEPERATE ELEMENT IN LIST


    ### UNCOMMENT IF ALL RAST AS STACK RASTER

    act_out = suppressWarnings(append(act_out, rast_points))

    ### UNCOMMENT IF ALL RAST AS STACK RASTER
  }

  if (normalize && norm_method == "range" && norm_group && length(data_iter) > 1){

    #max_val = max(terra::minmax(act_out)[2,], na.rm = TRUE)

    act_out = normalization(act_out, method = "range")

    # act_out = sapply(act_out, function(x) {
    #   if (!(all(is.na(terra::values(x))))) {
    #     terra::values(x) = BBmisc::normalize(terra::values(x), method = "range",
    #                                          margin = 2L,
    #                                          range = c(0, terra::minmax(x)[2] / max_val))
    #   }
    #   return(x)
    # })
  }


  if (!missing(env_data)){
    # project env_data to grid
    env_data_proj = terra::project(env_data, grid_rast)

    rast_env_points = act_out * env_data_proj
    output = rast_env_points
  } else {
    output = act_out
  }


  if (!missing(filepath)) { # save raster
      terra::writeRaster(output, filename = filepath)
      message(paste0("Saving output to ", filepath))
  }

  return(output)
}


