# library(twiGPS)
test_exposure_LS = function(data, x, y, NA_val, time_data, time_unit = "mins",
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
      data_proj = terra::na.omit(data_proj, time_cname)
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

  buff_list = list()
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

      #1st type of normalize - normalize time_elapsed without spatial reference
      if(normalize && (!norm_group || norm_method != "range" || length(data_iter) == 1)) {

        if (norm_method != "range" && norm_group) {
          message(paste0('Norm_method is "', norm_method, '" - norm_group is TRUE is applicable only for norm_method "range". Norm group argument ignored. Normalizing each group seperately'))
        }


        duration_line_id$time_elapsed = normalization(duration_line_id$time_elapsed,
                                                      method = norm_method)

      }
      #might use tidyterra for optimisation
      line_id_df = as.data.frame(traj_buff)
      buff_df = dplyr::left_join(line_id_df, duration_line_id, by = 'line_id')

      traj_buff$act = buff_df$time_elapsed

    }

    buff_list[length(buff_list) + 1] = list(traj_buff)

  }

  if(normalize && !missing(time_data) && norm_group && length(data_iter) > 1 &&
     norm_method == "range"){

    if (length(data_iter) > 1){

      max_buff = (do.call(rbind, buff_list))$act |> max(na.rm = TRUE)
      buff_list = sapply(buff_list, function(x) {
        if (!(all(is.na(x$act)))) {
          x$act = x$act / max_buff
        }
        return(x)
      })

    }
  }

  for (buff in buff_list) {

    if (!missing(env_data)){ #all of computation necessary only when calculating env_data
      # extract values and weights (area overlap) from grid for each cell of buffer
      exac_extr = exactextractr::exact_extract(grid_rast, sf::st_as_sf(buff), progress = FALSE)

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
      weight_buff = terra::merge(buff, traj_extract_line_id)

      weight_buff$end_weights = weight_buff$end_weights * weight_buff$act

      # rasterize results
      rast_segment = terra::rasterize(weight_buff, grid_rast, field = "end_weights", fun = "sum")
    } else {
      rast_segment = terra::rasterize(buff, grid_rast, field = "act", fun = "sum")
    }
    output = suppressWarnings(append(output, (rast_segment)))
  }

  if (!missing(filepath)) { # save raster
    terra::writeRaster(output, filename = filepath)
    message(paste0("Saving output to ", filepath))
  }


  return(output)

}
testthat::test_that("exposure_LS normalize range", {
  LS_test = exposure_LS(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "range",
                        time_data = dateTime, bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  LS =  test_exposure_LS(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "range",
                         time_data = dateTime, bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(LS) = "activity_space"
  testthat::expect_equal(LS_test,LS)
})

testthat::test_that("exposure_LS normalize center", {
  LS_test =  exposure_LS(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "center",
                                time_data = dateTime, bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  LS =  test_exposure_LS(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "center",
                         time_data = dateTime, bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(LS) = "activity_space"
  testthat::expect_equal(LS_test,LS)
})

testthat::test_that("exposure_LS normalize scale", {
  LS_test =  exposure_LS(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "scale",
                                time_data = dateTime, bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  LS =  test_exposure_LS(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "scale",
                         time_data = dateTime, bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(LS) = "activity_space"
  testthat::expect_equal(LS_test,LS)
})

testthat::test_that("exposure_LS normalize standardize", {
  LS_test =  exposure_LS(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "standardize",
                                time_data = dateTime, bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  LS =  test_exposure_LS(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "standardize",
                         time_data = dateTime, bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(LS) = "activity_space"
  testthat::expect_equal(LS_test,LS)
})

testthat::test_that("exposure_LS normalize range groups", {
  LS_test = exposure_LS(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "range",
                        time_data = dateTime, bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  LS =  test_exposure_LS(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "range",
                         time_data = dateTime, bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(LS) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
                "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
  testthat::expect_equal(LS_test,LS)
})

testthat::test_that("exposure_LS normalize range norm groups", {
  LS_test = exposure_LS(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "range",
                        time_data = dateTime, bandwidth = 200, norm_group = TRUE, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  LS =  test_exposure_LS(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "range",
                         time_data = dateTime, bandwidth = 200, norm_group = TRUE, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(LS) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
                "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
  testthat::expect_equal(LS_test,LS)
})

testthat::test_that("exposure_LS normalize center groups", {
  LS_test =  exposure_LS(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "center",
                                time_data = dateTime, bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  LS =  test_exposure_LS(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "center",
                         time_data = dateTime, bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(LS) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
                "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
  testthat::expect_equal(LS_test,LS)
})

testthat::test_that("exposure_LS normalize scale groups", {
  LS_test =  exposure_LS(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "scale",
                                time_data = dateTime, bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  LS =  test_exposure_LS(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "scale",
                         time_data = dateTime, bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(LS) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
                "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
  testthat::expect_equal(LS_test,LS)
})

testthat::test_that("exposure_LS normalize standardize groups", {
  LS_test =  exposure_LS(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "standardize",
                                time_data = dateTime, bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  LS =  test_exposure_LS(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "standardize",
                         time_data = dateTime, bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(LS) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
                "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
  testthat::expect_equal(LS_test,LS)
})

testthat::test_that("exposure_LS env_data", {
  ndvi_data = terra::rast(system.file("extdata/landsat_ndvi.tif", package = "twiGPS"))

  LS_test =  exposure_LS(data = geolife_sandiego ,coords = c("lon", "lat"), cellsize = 50,
                                time_data = dateTime, bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611", env_data = ndvi_data)
  LS =  test_exposure_LS(data = geolife_sandiego, x = lon, y = lat, cellsize = 50,
                         time_data = dateTime, bandwidth = 200, env_data = ndvi_data, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(LS) = "env_exposure"
  testthat::expect_equal(LS_test, LS)
})

