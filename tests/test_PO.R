library(twiGPS)

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

start_processing = function(data, x, y, NA_val, env_data, grid_extent,
                            input_crs, output_crs){

  # get spatial data
  if (all(class(data) == "data.frame")) {
    if (missing(input_crs)){ # implement so missing works
      input_crs = ""
    }
    if (!missing(NA_val)) {
      if (is.numeric(NA_val)){
        n_row = nrow(data)
        data = data[data[x] != NA_val & data[y] != NA_val,]
        warning(paste0("Removing rows containing default NA value ",
                       NA_val, " in x and y columns - removing ", n_row - nrow(data), " rows"))
      } else {
        warning("Invalid NA_val argument - NA_val is not numeric. Ignoring NA_val argument")
      }

    }
    if (any(is.na(data[,c(x,y)]))){
      n_row_NA = nrow(data)
      data = tidyr::drop_na(data, x, y)
      warning(paste0("Removing rows containing NA in x and y columns - removing", n_row_NA - nrow(data), " rows"))

    }
    data_points = terra::vect(x = data, geom = c(x, y), crs = input_crs)

  } else if (any(class(data) == "sf") && all(sf::st_geometry_type(data) == "POINT")){
    data_points = terra::vect(data)

    if (!missing(input_crs)){
      warning("Ignoring input_crs argument")
    }

  } else if (any(class(data) == "SpatVector") && terra::geomtype(data) == "points") {
    data_points = data
    if (!missing(input_crs)){
      warning("Ignoring input_crs argument")
    }
  } else {
    stop("Object data is neither data.frame, sf nor SpatVector class")
  }

  data_crs = terra::crs(data_points)

  # crs
  if (!missing(output_crs)){
    if (data_crs == "") {
      message("Setting data crs to output_crs")
      terra::crs(data_points) = output_crs
      data_proj = data_points
    } else {
      message("Projecting data to output_crs")
      data_proj = terra::project(data_points, output_crs)
    }
  } else if (!missing(grid_extent) && class(grid_extent) == "SpatRaster") {
    if (data_crs == "") {
      message("Setting data crs to grid_extent crs")
      terra::crs(data_points) = terra::crs(grid_extent)
      data_proj = data_points
    } else {
      message("Projecting data to grid_extent crs")
      data_proj = terra::project(data_points, grid_extent)
    }
  } else if (data_crs != "") { # any invalid/empty crs
    data_proj = data_points
  } else if (!missing(env_data)) {
    message("Setting data crs to environmental data crs")
    terra::crs(data_points) = terra::crs(env_data)
    data_proj = data_points
  } else {
    message("No crs specified")
    data_proj = data_points
  }
  # crop to grid_extent
  if (!missing(grid_extent)){
    if (class(grid_extent) == "SpatExtent"){
      data_proj = terra::crop(data_proj, grid_extent)
    } else {
      if (terra::crs(data_proj) != terra::crs(grid_extent)){
        grid_rast = terra::project(grid_rast, terra::crs(data_proj))
      }
      data_proj = terra::crop(data_proj, terra::ext(grid_extent))
    }
  }

  return(data_proj)

}
calc_grid = function(x, bandwidth, cellsize, env_data, grid_extent, is_LS = FALSE){

  if (!missing(grid_extent)){ # grid_extent as ext of grid rast

    extent = grid_extent
  } else { # calculate extent
    extent = terra::ext(x)
    if (!missing(bandwidth)) { # for KDE/DR expand extent by bandwidth
      if(!is_LS){
        extent = c(terra::xmin(extent) - bandwidth, terra::xmax(extent) + bandwidth,
                   terra::ymin(extent) - bandwidth, terra::ymax(extent) + bandwidth)
      } else { # for LS in lon/lat crs bandwidth is still in meters (terra::buffer)
        temp_buff = terra::buffer(x, bandwidth)
        extent = terra::ext(temp_buff)
      }

    }
  }


  if(!missing(cellsize) && length(cellsize) == 1 && is.numeric(cellsize) && cellsize > 0) { # cellsize included

    if (suppressWarnings(!is.na(terra::is.lonlat(x)) && terra::is.lonlat(x))){
      # crs units in degrees
      # is.na if empty crs to skip error
      warning("Cellsize is not stable - cells are not rectangular")

      if (!missing(bandwidth) && bandwidth > 0.1 && !is_LS) {
        message(paste0("CRS is in lontitude/latitude and bandwidth is ", bandwidth,
                       " - bandwidth is calculated in CRS units. Is bandwidth in correct unit?"))
      }

      dist_lon = geosphere::distm(c(extent[1], extent[3]), c(extent[2], extent[3]),
                                  fun = geosphere::distHaversine)
      dist_lat = geosphere::distm(c(extent[1], extent[3]), c(extent[1], extent[4]),
                                  fun = geosphere::distHaversine)

      # number of cells
      x_cells = (dist_lon / cellsize) |> as.integer()
      y_cells = (dist_lat / cellsize) |> as.integer()

      # empty_rast for units in degrees
      grid_rast = terra::rast(crs = terra::crs(x), nrows=y_cells,
                              ncols=x_cells, extent = extent)

    } else {
      grid_rast = terra::rast(crs = terra::crs(x), extent = extent,
                              resolution = cellsize)
    }

  } else if  (!missing(env_data) && any(class(env_data) == "SpatRaster")){ #if incorrect cellsize and env_data exists
    warning("Cellsize invalid or not set - using cellsize from env_data")

    if (terra::linearUnits(env_data) == terra::linearUnits(x)){
      # project env data with the same cellsize
      env_data_proj = terra::project(env_data, terra::crs(x), res = terra::res(env_data)[1])
      # if env_data in lat/lon only quadratic cells

    } else {
      env_data_proj = terra::project(env_data, terra::crs(x))
    }
    if (missing(grid_extent)){
      # extent modified to preserve cellsize
      grid_rast = terra::crop(env_data_proj, extent)
    } else {
      # if grid_extent specified and missing cellsize grid_extent is prior to
      # cellsize and cellsize is slightly modified
      grid_rast = env_data_proj
      terra::ext(grid_rast) = extent
    }
  } else {
    stop("Invalid or missing cellsize - provide valid cellsize or env_data")
  }
  return(grid_rast)
}
env_vect = function(env, env_buff, env_field, grid){

  if (!missing(env_buff)) { # create buffer around vector data
    if (length(env_buff) == 1 && is.numeric(env_buff) && env_buff > 0){
      env_proj = terra::project(env, terra::crs(grid))
      ## TERRA::BUFFER SOMETIMES CRUSHES WITH BIG LINE VECTOR DATASETS
      #env = terra::buffer(env_proj, env_buff)
      ## TEMPORARY FIX: USE SF::ST_BUFFER
      env = env_proj |> sf::st_as_sf() |> sf::st_buffer(env_buff) |> terra::vect()

    } else {
      warning("Invalid env_buff argument - ignoring creation of buffer")
    }
  }

  if (!missing(env_field)){ # choose field to rasterize
    env_f_enq = rlang::enquo(env_field)
    if (rlang::quo_name(env_f_enq) %in% terra::names(env)){
      env = terra::rasterize(env, grid, field = rlang::quo_name(env_f_enq), fun = "sum")
    } else {
      warning("Invalid env_field argument - env_field is not a column name in data. Ignoring env_field argument")
      env = terra::rasterize(env, grid, fun = "sum")
    }
  } else {
    env = terra::rasterize(env, grid, fun = "sum")
  }

  return(env)
}
normalization = function(data, method, range = c(0, 1)){
  if (inherits(data, "SpatRaster") && diff(c(min(terra::minmax(data)[1,], na.rm = TRUE),
                                             max(terra::minmax(data)[2,], na.rm = TRUE))) == 0 ) {
    switch(method,
           center = terra::scale(data, center = TRUE, scale = FALSE),
           # range = (data - terra::minmax(data, na.rm = TRUE)[1,]) /
           #   diff(terra::minmax(data, na.rm = TRUE)) * diff(range) + range[1L],
           #without range arg
           range = data / max(terra::minmax(data)[2,], na.rm = TRUE),
           standardize = terra::scale(data, center = TRUE, scale = FALSE),
           scale = data)
  } else if (inherits(data, "numeric") && length(unique(data[!is.na(data)])) == 1){
    switch(method,
           center = terra::scale(data, center = TRUE, scale = FALSE),
           # range = (data - terra::minmax(data, na.rm = TRUE)[1,]) /
           #   diff(terra::minmax(data, na.rm = TRUE)) * diff(range) + range[1L],
           #without range arg
           range = data / max(data, na.rm = TRUE),
           standardize = terra::scale(data, center = TRUE, scale = FALSE),
           scale = data)
  } else if (inherits(data, "numeric") && length(unique(data[!is.na(data)])) != 1) {
    switch(method,
           # range = (data - terra::minmax(data, na.rm = TRUE)[1,]) /
           #   diff(terra::minmax(data, na.rm = TRUE)) * diff(range) + range[1L],
           #without range arg
           range = data / max(data, na.rm = TRUE),
           standardize = scale(data, center = TRUE, scale = TRUE),
           center = scale(data, center = TRUE, scale = FALSE),
           scale = scale(data, center = FALSE, scale = stats::sd(data, na.rm = TRUE))
    )
  } else {
    switch(method,
           # range = (data - terra::minmax(data, na.rm = TRUE)[1,]) /
           #   diff(terra::minmax(data, na.rm = TRUE)) * diff(range) + range[1L],
           #without range arg
           range = data / max(terra::minmax(data)[2,], na.rm = TRUE),
           standardize = terra::scale(data, center = TRUE, scale = TRUE),
           center = terra::scale(data, center = TRUE, scale = FALSE),
           scale = terra::scale(data, center = FALSE, scale = terra::global(data, "sd", na.rm = TRUE)[[1]])
    )
  }
}
testthat::test_that("exposure_PO normalize range", {
  PO_test = exposure_PO(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "range",
                   input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  PO =  test_exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "range",
                    input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  testthat::expect_equal(PO_test, PO)
})

testthat::test_that("exposure_PO normalize center", {
  PO_test = exposure_PO(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "center",
                        input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  PO =  test_exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "center",
                         input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  testthat::expect_equal(PO_test, PO)
})

testthat::test_that("exposure_PO normalize scale", {
  PO_test = exposure_PO(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "scale",
                        input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  PO =  test_exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "scale",
                         input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  testthat::expect_equal(PO_test, PO)
})

testthat::test_that("exposure_PO normalize standardize", {
  PO_test = exposure_PO(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "standardize",
                        input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  PO =  test_exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "standardize",
                         input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  testthat::expect_equal(PO_test, PO)
})


testthat::test_that("exposure_PO normalize range groups", {
  PO_test = exposure_PO(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "range",
                        group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  PO =  test_exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "range",
                         group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  testthat::expect_equal(PO_test, PO)
})

testthat::test_that("exposure_PO normalize range groups", {
  PO_test = exposure_PO(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "range",
                        norm_group = TRUE, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  PO =  test_exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "range",
                         norm_group = TRUE, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  testthat::expect_equal(PO_test, PO)
})

testthat::test_that("exposure_PO normalize center groups", {
  PO_test = exposure_PO(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "center",
                                group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  PO =  test_exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "center",
                         group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  testthat::expect_equal(PO_test, PO)
})

testthat::test_that("exposure_PO normalize scale groups", {
  PO_test = exposure_PO(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "scale",
                                group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  PO =  test_exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "scale",
                         group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  testthat::expect_equal(PO_test, PO)
})

testthat::test_that("exposure_PO normalize standardize groups", {
  PO_test = exposure_PO(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "standardize",
                                group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  PO =  test_exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "standardize",
                         group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  testthat::expect_equal(PO_test, PO)
})

testthat::test_that("exposure_PO env_data", {
  ndvi_data = terra::rast(system.file("extdata/landsat_ndvi.tif", package = "twiGPS"))

  PO_test = twiGPS::exposure_PO(data = twiGPS::geolife_sandiego, coords = c("lon", "lat"), cellsize = 50,
                                 input_crs = "EPSG:4326", output_crs = "EPSG:32611", env_data = ndvi_data)
  PO =  test_exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50,
                         env_data = ndvi_data, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  testthat::expect_equal(PO_test, PO)
})
