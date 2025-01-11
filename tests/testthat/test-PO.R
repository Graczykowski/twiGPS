test_exposure_PO = function(data, x, y, NA_val, cellsize, group_split, env_data, env_field, env_buff,
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
  output = list()


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
      if (norm_method == "range"){
        rast_points = terra::subst(rast_points, from = NA, to = 0)
        rast_points = normalization(rast_points, method = norm_method)
        rast_points = terra::subst(rast_points, from = 0, to = NA)
      } else {
        rast_points = normalization(rast_points, method = norm_method)
      }
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
    # if (length(data_iter) > 1){ # update filepath
    #
    #   mapply(function(x, y) {
    #     terra::writeRaster(x, filename = y)
    #     message(paste0("Saving output to ", y))
    #   }, x = output, y = file_vect)
    #
    # } else {

    terra::writeRaster(end_rast, filename = filepath)
    message(paste0("Saving output to ", filepath))
    # }
  }

  return(output)
}

testthat::test_that("exposure_PO normalize range", {
  PO_test = exposure_PO(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "range",
                   input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  PO =  test_exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "range",
                    input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(PO) = "activity_space"
  testthat::expect_equal(PO_test, PO)
})

testthat::test_that("exposure_PO normalize center", {
  PO_test = exposure_PO(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "center",
                        input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  PO =  test_exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "center",
                         input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(PO) = "activity_space"
  testthat::expect_equal(PO_test, PO)
})

testthat::test_that("exposure_PO normalize scale", {
  PO_test = exposure_PO(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "scale",
                        input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  PO =  test_exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "scale",
                         input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(PO) = "activity_space"
  testthat::expect_equal(PO_test, PO)
})

testthat::test_that("exposure_PO normalize standardize", {
  PO_test = exposure_PO(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "standardize",
                        input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  PO =  test_exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "standardize",
                         input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(PO) = "activity_space"
  testthat::expect_equal(PO_test, PO)
})


testthat::test_that("exposure_PO normalize range groups", {
  PO_test = exposure_PO(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "range",
                        group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  PO =  test_exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "range",
                         group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(PO) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
                "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
  testthat::expect_equal(PO_test, PO)
})

testthat::test_that("exposure_PO normalize range groups", {
  PO_test = exposure_PO(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "range",
                        norm_group = TRUE, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  PO =  test_exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "range",
                         norm_group = TRUE, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(PO) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
                "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
  testthat::expect_equal(PO_test, PO)
})

testthat::test_that("exposure_PO normalize center groups", {
  PO_test = exposure_PO(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "center",
                                group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  PO =  test_exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "center",
                         group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(PO) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
                "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
  testthat::expect_equal(PO_test, PO)
})

testthat::test_that("exposure_PO normalize scale groups", {
  PO_test = exposure_PO(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "scale",
                                group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  PO =  test_exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "scale",
                         group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(PO) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
                "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
  testthat::expect_equal(PO_test, PO)
})

testthat::test_that("exposure_PO normalize standardize groups", {
  PO_test = exposure_PO(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "standardize",
                                group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  PO =  test_exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "standardize",
                         group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(PO) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
                "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
  testthat::expect_equal(PO_test, PO)
})

testthat::test_that("exposure_PO env_data", {
  ndvi_data = terra::rast(system.file("extdata/landsat_ndvi.tif", package = "twiGPS"))

  PO_test = twiGPS::exposure_PO(data = twiGPS::geolife_sandiego, coords = c("lon", "lat"), cellsize = 50,
                                 input_crs = "EPSG:4326", output_crs = "EPSG:32611", env_data = ndvi_data)
  PO =  test_exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50,
                         env_data = ndvi_data, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(PO) = "env_exposure"
  testthat::expect_equal(PO_test, PO)
})
