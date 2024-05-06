#' Density Ranking exposure
#'
#' @description
#' Density Ranking method activity space and environmental exposure.
#' Using spat_dr() function based on Silverman, B. W. 1986 (KDE) and DR function from Yen-Chi Chen's repository yenchic/density_ranking.
#' In order to receive activity space ignore env_data argument.
#'
#' @param data Data.frame, SpatVector points or sf data.frame containing only POINTS.
#' @param x Data masking. x or longitude coordinates column name if data is a data.frame.
#' @param y Data masking. y or latitude coordinates column name if data is a data.frame.
#' @param NA_val Numeric. Value in x and y marked as NA if data is a data.frame.
#' @param cellsize Positive numeric. Size of raster cells in meters if output has a longtitude/latitude CRS or in the units of output CRS.
#' @param group_split Data masking. Column of data based on which it is grouped and split. For n groups output will be n SpatRasters.
#' @param bandwidth Positive numeric. Bandwidth in units of output Coordinate Reference System.
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
#'
#' @return SpatRaster
#'
#' @examples
#'
#' # activity space
#' exposure_DR(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, bandwidth = 200,
#'  normalize = TRUE, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
#'
#' # split by date
#' exposure_DR(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, group_split = date,
#'    bandwidth = 200, normalize = TRUE, norm_group = FALSE, input_crs = "EPSG:4326",
#'    output_crs = "EPSG:32611")
#'
#' # split by date and define extent
#' extent = c(478000, 484000, 3618000, 3627000)
#'
#' exposure_DR(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, group_split = date,
#'   bandwidth = 200, normalize = TRUE,norm_group = TRUE, grid_extent = extent,
#'   input_crs = "EPSG:4326", output_crs = "EPSG:32611")
#'
#' #environmental exposure
#'
#' ndvi_data = terra::rast(system.file("extdata/landsat_ndvi.tif", package = "twsagps"))
#'
#' exposure_DR(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, bandwidth = 200,
#'   env_data = ndvi_data, normalize = FALSE, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
#'
#' # environmental exposure - use rast grid and split by date
#' exposure_DR(data = geolife_sandiego, x = lon, y = lat, group_split = date,
#'   bandwidth = 200, env_data = ndvi_data, normalize = TRUE, norm_group = TRUE,
#'   grid_extent = ndvi_data, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
#'
#' @export


exposure_DR = function(data, x, y, NA_val, cellsize, group_split, bandwidth, env_data,
                       env_field, env_buff, normalize = FALSE, norm_method = "range",
                       norm_group = FALSE, grid_extent, input_crs, output_crs, filepath){


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

  #handle bandwidth
  if (missing(bandwidth)){
    stop("Missing bandwidth argument. Provide valid bandwidth")
  } else if (length(bandwidth) != 1 || !is.numeric(bandwidth) || bandwidth <= 0){
    stop("Invalid bandwidth argument - bandwidth is neither positive nor single numeric value")
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

  # handle normalize and norm_group
  if (is.na(as.logical(normalize))){
    stop("Invalid normalize argument - normalize argument cannot be interpreted as boolean")
  } else {
    normalize = as.logical(normalize)
  }

  if (is.na(as.logical(norm_group))){
    stop("Invalid norm_group argument - norm_group argument cannot be interpreted as boolean")
  } else {
    norm_group = as.logical(norm_group)
  }

  # handle norm_method
  if (!norm_method %in% c("center", "scale", "standardize", "range") && normalize == TRUE) {
    warning('Invalid norm_method - applying default normalization method "range"')
    norm_method = "range"
  }

  # get spatial data with correct crs
  if (!missing(data)){
    if (all(class(data) == "data.frame")){
      if (all(!c(missing(x), missing(y)))){
        x_enq = rlang::enquo(x)
        y_enq = rlang::enquo(y)

        if (all(c(rlang::quo_name(y_enq), rlang::quo_name(y_enq)) %in% colnames(data))){

          data_proj = start_processing(data = data, x = rlang::quo_name(x_enq),
                                       y = rlang::quo_name(y_enq), NA_val = NA_val,
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

  if (!missing(grid_extent) && any(class(grid_extent) == "SpatRaster")){
    if (terra::crs(data_proj) != terra::crs(grid_extent)){
      grid_extent = terra::project(grid_extent, terra::crs(data_proj))
    }
    grid_rast = grid_extent
  } else {
    grid_rast = calc_grid(x = data_proj, cellsize = cellsize, env_data = env_data,
                          grid_extent = grid_extent, bandwidth = bandwidth)
  }

  # if env_data is vector data - create optional buffer and rasterize to grid raster
  if (!missing(env_data) && any(class(env_data) == "SpatVector")) {

    env_data = env_vect(env = env_data, env_buff = env_buff,
                        env_field = env_field, grid = grid_rast)
  }


  if (missing(group_split)) {
    data_iter = list(data_proj) # only one item for for loop
  } else {
    enq_group_split = rlang::enquo(group_split)
    if (rlang::quo_name(enq_group_split) %in% terra::names(data_proj)){
      data_iter = terra::split(data_proj, rlang::quo_name(enq_group_split)) # split data_proj by group_split
      message(paste0("Data split by group into ", length(data_iter), " items"))
    } else {
      data_iter = list(data_proj)
      warning("Invalid group_split argument - group_split is not a column in data. Data not split")
    }

  }

  if (length(data_iter) > 1 && !missing(filepath)) {
    # list output

    # for iterating file names
    file_ext = stringi::stri_extract(filepath, regex = "\\.(\\w+)$")
    file_no_ext = substr(filepath, 1, nchar(filepath) - nchar(file_ext))

    file_vect = rep(file_no_ext, length.out = length(data_iter))
    file_vect = paste0(file_vect, "_", seq_along(file_vect), file_ext)
  }

  act_out = list()
  output = list()

  for (data_i in data_iter){
    # if each group should have seperate extent then output is a list rasts
    # if all groups should have same extent then output is rast with n layers

    ### UNCOMMENT IF EVERY RAST SHOULD HAVE SEPERATE EXTENT

    if (missing(grid_extent)) {
      #new ext for each group
      # get extent
      group_extent = terra::ext(data_i)
      # new extent - expanded extent by bandwidth
      new_group_extent = c(terra::xmin(group_extent) - bandwidth,
                           terra::xmax(group_extent) + bandwidth,
                           terra::ymin(group_extent) - bandwidth,
                           terra::ymax(group_extent) + bandwidth)

      # crop ext of each rast
      grid_crop = terra::crop(grid_rast, new_group_extent)

      dr_rast = spat_dr(data_i, grid_crop, bandwidth)
    } else {
      dr_rast =  spat_dr(data_i, grid_rast, bandwidth)
    }
    ### UNCOMMENT IF EVERY RAST SHOULD HAVE SEPERATE EXTENT


    ### UNCOMMENT IF EVERY RAST SHOULD HAVE SAME EXTENT

    #dr_rast = spat_dr(data_i, grid_rast, bandwidth)

    ### UNCOMMENT IF EVERY RAST SHOULD HAVE SAME EXTENT

    if (normalize == TRUE && (norm_group == FALSE || norm_method != "range" || length(data_iter) == 1)){
      if (norm_method != "range" && norm_group == TRUE) {
        message(paste0('Norm_method is "', norm_method, '" - norm_group is TRUE is applicable only for norm_method "range". Norm group argument ignored. Normalizing each group seperately'))
      }
      # calculate normalization to 0-1 range
      terra::values(dr_rast) = BBmisc::normalize(terra::values(dr_rast),
                                                 method = norm_method, margin = 2L)

    }

    dr_rast = terra::subst(dr_rast, from = 0, to = NA)

    ### UNCOMMENT IF EVERY RAST SHOULD BE A SEPERATE ELEMENT IN LIST

    act_out[length(act_out) + 1] = as.list(dr_rast)

    ### UNCOMMENT IF EVERY RAST SHOULD BE A SEPERATE ELEMENT IN LIST


    ### UNCOMMENT IF ALL RAST AS STACK RASTER

    #act_out = append(act_out, kde_rast)

    ### UNCOMMENT IF ALL RAST AS STACK RASTER
  }

  if (normalize == TRUE && norm_method == "range" && norm_group == TRUE && length(data_iter) > 1){

    max_val = max(sapply(act_out, terra::minmax)[2,], na.rm = TRUE)

    act_out = sapply(act_out, function(x) {
      if (!(all(is.na(terra::values(x))))) {
      terra::values(x) = BBmisc::normalize(terra::values(x), method = "range",
                                           margin = 2L,
                                           range = c(0, terra::minmax(x)[2] / max_val))
      }
      return(x)
    })
  }

  for (out in act_out) {
    if (!missing(env_data)){ # calculate exposure
      # project env_data to grid
      env_data_proj = terra::project(env_data, out)
      # if same ext for groups out -> grid_rast and this code chunk out of loop

      #env_data_resamp = terra::resample(env_data_proj, out)

      rast_env_dr = out * env_data_proj
      end_rast = rast_env_dr
    } else {
      end_rast = out
    }
    if (length(data_iter) == 1) { # if one group or no group output is not a list
      output = end_rast
    } else {
      output[length(output) + 1] = as.list(end_rast)
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
