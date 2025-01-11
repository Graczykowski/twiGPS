# library(twiGPS)

test_exposure_KDE = function(data, x, y, NA_val, cellsize, group_split, bandwidth, env_data,
                        env_field, env_buff, normalize = FALSE, norm_method = "range",
                        norm_group = FALSE, grid_extent, input_crs, output_crs, filepath){

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
    #   #new ext for each group
    #   # get extent
    #   group_extent = terra::ext(data_i)
    #   # new extent - expanded extent by bandwidth
    #   new_group_extent = c(terra::xmin(group_extent) - bandwidth,
    #                        terra::xmax(group_extent) + bandwidth,
    #                        terra::ymin(group_extent) - bandwidth,
    #                        terra::ymax(group_extent) + bandwidth)
    #
    #   # crop ext of each rast
    #   grid_crop = terra::crop(grid_rast, new_group_extent)
    #
    #   kde_rast = spat_kde(data_i, grid_crop, bandwidth)
    # } else {
    #   kde_rast =  spat_kde(data_i, grid_rast, bandwidth)
    # }

    ### UNCOMMENT IF EVERY RAST SHOULD HAVE SEPERATE EXTENT


    ### UNCOMMENT IF EVERY RAST SHOULD HAVE SAME EXTENT

    kde_rast = spat_kde(data_i, grid_rast, bandwidth)

    ### UNCOMMENT IF EVERY RAST SHOULD HAVE SAME EXTENT

    if (normalize && (!norm_group || norm_method != "range" || length(data_iter) == 1)){
      if (norm_method != "range" && norm_group) {
        message(paste0('Norm_method is "', norm_method, '" - norm_group is TRUE is applicable only for norm_method "range". Norm group argument ignored. Normalizing each group seperately'))
      }
      # calculate normalization
      kde_rast = normalization(kde_rast, method = norm_method)

    }

    kde_rast = terra::subst(kde_rast, from = 0, to = NA)


    ### UNCOMMENT IF EVERY RAST SHOULD BE A SEPERATE ELEMENT IN LIST

    # act_out[length(act_out) + 1] = as.list(kde_rast)

    ### UNCOMMENT IF EVERY RAST SHOULD BE A SEPERATE ELEMENT IN LIST


    ### UNCOMMENT IF ALL RAST AS STACK RASTER

    act_out = append(act_out, kde_rast)

    ### UNCOMMENT IF ALL RAST AS STACK RASTER
  }

  if (normalize && norm_method == "range" && norm_group && length(data_iter) > 1){

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
    terra::writeRaster(output, filename = filepath)
    message(paste0("Saving output to ", filepath))
  }

  return(output)

}


testthat::test_that("exposure_KDE normalize range", {
  KDE_test =  exposure_KDE(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "range",
                        bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  KDE =  test_exposure_KDE(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "range",
                         bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(KDE) = "activity_space"
  testthat::expect_equal(KDE_test,KDE)
})

testthat::test_that("exposure_KDE normalize center", {
  KDE_test =  exposure_KDE(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "center",
                                bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  KDE =  test_exposure_KDE(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "center",
                         bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(KDE) = "activity_space"
  testthat::expect_equal(KDE_test,KDE)
})

testthat::test_that("exposure_KDE normalize scale", {
  KDE_test =  exposure_KDE(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "scale",
                                bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  KDE =  test_exposure_KDE(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "scale",
                         bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(KDE) = "activity_space"
  testthat::expect_equal(KDE_test,KDE)
})

testthat::test_that("exposure_KDE normalize standardize", {
  KDE_test =  exposure_KDE(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "standardize",
                                bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  KDE =  test_exposure_KDE(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "standardize",
                         bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(KDE) = "activity_space"
  testthat::expect_equal(KDE_test,KDE)
})


testthat::test_that("exposure_KDE normalize range groups", {
  KDE_test = exposure_KDE(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "range",
                        bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  KDE =  test_exposure_KDE(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "range",
                         bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(KDE) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
                "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
  testthat::expect_equal(KDE_test,KDE)
})

testthat::test_that("exposure_KDE normalize range groups", {
  KDE_test = exposure_KDE(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "range",
                        bandwidth = 200, norm_group = TRUE, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  KDE =  test_exposure_KDE(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "range",
                         bandwidth = 200, norm_group = TRUE, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(KDE) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
                 "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
  testthat::expect_equal(KDE_test,KDE)
})

testthat::test_that("exposure_KDE normalize center groups", {
  KDE_test =  exposure_KDE(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "center",
                                bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  KDE =  test_exposure_KDE(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "center",
                         bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(KDE) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
                 "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
  testthat::expect_equal(KDE_test,KDE)
})

testthat::test_that("exposure_KDE normalize scale groups", {
  KDE_test =  exposure_KDE(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "scale",
                                bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  KDE =  test_exposure_KDE(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "scale",
                         bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(KDE) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
                 "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
  testthat::expect_equal(KDE_test,KDE)
})

testthat::test_that("exposure_KDE normalize standardize groups", {
  KDE_test =  exposure_KDE(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "standardize",
                                bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  KDE =  test_exposure_KDE(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "standardize",
                         bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(KDE) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
                 "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
  testthat::expect_equal(KDE_test,KDE)
})

testthat::test_that("exposure_KDE env_data", {
  ndvi_data = terra::rast(system.file("extdata/landsat_ndvi.tif", package = "twiGPS"))

  KDE_test =  exposure_KDE(data = geolife_sandiego ,coords = c("lon", "lat"), cellsize = 50,
                                bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611", env_data = ndvi_data)
  KDE =  test_exposure_KDE(data = geolife_sandiego, x = lon, y = lat, cellsize = 50,
                         bandwidth = 200, env_data = ndvi_data, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(KDE) = "env_exposure"

  testthat::expect_equal(KDE_test, KDE)
})
