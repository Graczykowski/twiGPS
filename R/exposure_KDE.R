#' Calculate Kernel Density Estimation exposure of point spatial data
#'
#' @description
#' `exposure_KDE()` calculates activity space and environmental exposure using Kernel Density Estimation method.
#' Kernel Density Estimation function is based on Silverman's quartic kernel:
#' \deqn{
#' \hat{f}(x) = \frac{1}{h^2} \sum^n_{i=1} K \{\frac{1}{h}(x - X_i)\}
#' }
#' and
#' \deqn{
#' K_2(x) = \left\{ \begin{array}{rcl} 3\pi^{-1}(1-x^Tx)^2 & if & x^Tx > 1 \\ 0 & \ & otherwise \end{array}\right.
#' }
#' where `h` is a bandwidth of kernel (`bandwidth` argument). Formula calculates output in intensity units.
#' Setting `env_data` argument determines whether environmental exposure or activity space is computed.
#'
#' For calculating statistics of output raster use [exposure_stats()]
#'
#' @param data Data.frame, SpatVector points or sf data.frame containing only POINTS.
#' @param coords Character. Vector of column names of x and y coordinates if `data` is a data.frame.
#' @param bandwidth Positive numeric. Bandwidth in units of output's CRS.
#' @param cellsize Positive numeric. Size of raster cells in meters if output has a longtitude/latitude CRS or in the units of output CRS.
#' @param env_data Stars, SpatRaster, SpatVector or sf. Spatial environmental data. Activity space is calculated when not set. When argument is a SpatVector or sf object vector data is rasterized to output raster using `"sum"` function.
#' @param output_crs Character or terra crs object. CRS of output. If not set and `env_data` is a SpatRaster `env_data`'s CRS is used.
#' @param input_crs Character or terra crs object. CRS of `data`'s coordinates if `data` is a data.frame.
#' @param grid_extent Stars, SpatRaster, SpatExtent, sf bbox object or numeric vector of 4 length `c(xmin, xmax, ymin, ymax)`. If stars or SpatRaster `grid_extent` is output's grid and `cellsize` argument is ignored. If SpatExtent, sf bbox object or vector `grid_extent` is output's extent. Then, if `cellsize` is set it is preserved at the cost of extent.
#' @param normalize Character. If set activity space SpatRaster is normalized with specified method. Four methods - "center", "scale", "standardize" and "range". See [BBmisc::normalize()].
#' @param group_split Data masking. Column of data based on which it is grouped and split. For n groups output will be SpatRasters with n layers.
#' @param norm_group Boolean. Applicable only When normalize is "range" and group_split is set. If FALSE each layer of SpatRaster is normalized seperately. If TRUE all layer are normalized to layer with highest max value (default FALSE).
#' @param filepath Character. Output filename. See [terra::writeRaster()].
#' @param NA_val Numeric. Value in x and y marked as NA if `data` is a data.frame.
#' @param env_field Data masking. Column of `env_data` that `env_data` will be rasterized on. Ignored if `env_data` not a SpatVector or sf class.
#' @param env_buff Positive numeric. Optional buffer around SpatVector/sf `env_data` in meters if output has a longtitude/latitude CRS or in the units of the CRS. Ignored if `env_data` not a SpatVector or sf class. See [terra::buffer()].
#' @param verbose Boolean. If FALSE amount of output is reduced (default TRUE).
#'
#' @return SpatRaster
#'
#' @references Jankowska, Marta M., Jiue-An Yang, Nana Luo, Chad Spoon, and Tarik Benmarhnia. 2023. “Accounting for Space, Time, and Behavior Using GPS Derived Dynamic Measures of Environmental Exposure.” Health & Place 79 (January): 102706. https://doi.org/10.1016/j.healthplace.2021.102706.
#'
#' Silverman, B. W. 1986. Density Estimation for Statistics and Data Analysis. Vol. 37. 1. Chapman Hall.
#' @examples
#' # activity space for data.frame 'data'
#' exposure_KDE(data = geolife_sandiego, coords = c("lon", "lat"),
#'   bandwidth = 200, cellsize = 50, output_crs = "EPSG:32611",
#'   input_crs = "EPSG:4326")
#'
#' # SpatVector data
#' geolife_vect = terra::vect(geolife_sandiego,
#'                            geom = c("lon", "lat"),
#'                            crs = "EPSG:4326") |>
#'                            terra::project("EPSG:32611")
#'
#' # normalize "range"
#' exposure_KDE(data = geolife_vect, bandwidth = 200,
#'   cellsize = 50, normalize = "range")
#'
#' # normalize "standardize"
#' exposure_KDE(data = geolife_vect, bandwidth = 200,
#'   cellsize = 50, normalize = "standardize")
#'
#' # split by date
#' exposure_KDE(data = geolife_vect, bandwidth = 200,
#'   cellsize = 50, normalize = "range", group_split = date)
#'
#' # split by date, normalize by group and define extent
#' extent = c(478000, 484000, 3618000, 3627000)
#'
#' exposure_KDE(data = geolife_vect, bandwidth = 200,
#'   cellsize = 50, grid_extent = extent, normalize = "range",
#'   group_split = date, norm_group = TRUE)
#'
#' # environmental exposure
#' ndvi_data = terra::rast(system.file("extdata/landsat_ndvi.tif",
#'                                      package = "twiGPS"))
#'
#' exposure_KDE(data = geolife_vect, bandwidth = 200,
#'  cellsize = 50, env_data = ndvi_data)
#'
#' # environmental exposure - use rast grid and split by date
#' exposure_KDE(data = geolife_vect, bandwidth = 200,
#'   env_data = ndvi_data, grid_extent = ndvi_data,
#'   normalize = "range", group_split = date)
#'
#' @seealso [exposure_PO()], [exposure_DR()], [exposure_LS()]
#' @importFrom rlang .data
#' @export
exposure_KDE = function(data, coords, bandwidth, cellsize, env_data, output_crs,
                        input_crs, grid_extent, normalize, group_split,
                        norm_group = FALSE, filepath, NA_val, env_field,
                        env_buff, verbose = TRUE){


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

  act_out = list()

  for (data_i in data_iter){

    kde_rast = spat_kde(data_i, grid_rast, bandwidth)

    if (!missing(group_split) && enq_group_split %in% names(data_proj)){
      r_name = unique(data_i[enq_group_split][[1]])
    } else {
      r_name = "activity_space"
    }
    names(kde_rast) = r_name

    if (!missing(normalize) && (!norm_group || normalize != "range" || length(data_iter) == 1)){
      if (normalize != "range" && norm_group && verbose) {
        message(paste0('Normalization method is "', normalize, '" - \'norm_group\' = TRUE is applicable only for normalization method "range". \'Norm_group\' argument ignored, each group is normalized seperately'))
      }
      # calculate normalization
      kde_rast = normalization(kde_rast, method = normalize)

    }

    kde_rast = terra::subst(kde_rast, from = 0, to = NA)

    act_out = suppressWarnings(append(act_out, kde_rast))
  }

  if (!missing(normalize) && normalize == "range" && norm_group && length(data_iter) > 1){

    act_out = normalization(act_out, method = "range")

  }

  if (!missing(env_data)){
    # project env_data to grid
    env_data_proj = terra::project(env_data, grid_rast)

    rast_env_kde = act_out * env_data_proj

    if (missing(group_split) || !enq_group_split %in% names(data_proj)){
      names(rast_env_kde) = "env_exposure"
    }

    output = rast_env_kde
  } else {
    output = act_out
  }



  if (!missing(filepath)) { # save raster
    terra::writeRaster(output, filename = filepath)
    if (verbose) {
      message(paste0("Saving output to ", filepath))
    }
  }

  return(output)

}
