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
trajectories_fun = function(data){

  t_name = c(names(data))[grepl('time', c(names(data)))]
  lst_name = paste(t_name[which.max(nchar(t_name))], "_1", sep = "")

  trajectories_out = data |>
    sf::st_as_sf() |> # change to sf object to change points to linestring later
    #filter to the study id defined by the argument
    #define start and end points of line
    dplyr::mutate(
      line_id = dplyr::row_number(),#an id for each "line segment"
      x_start= sf::st_coordinates(geometry)[,1],
      y_start= sf::st_coordinates(geometry)[,2],
      x_end = dplyr::lead(x_start),
      y_end = dplyr::lead(y_start)
    ) |>
    dplyr::ungroup() |>
    sf::st_set_geometry(NULL) |>
    #exclude the last observation, which has no "lead", and will be missing.
    dplyr::filter(is.na(x_end)==FALSE) |>
    # Make the data long form so that each point has two observations
    tidyr::pivot_longer(
      # select variables to pivot longer.
      cols = c(x_start, y_start, x_end, y_end),
      #value goes to "x/y", and time goes to "_start/end"
      names_to = c(".value", lst_name),
      names_repair = "unique",
      names_sep = "_"#the separator for the column name
    ) |>
    # create sf object once again
    sf::st_as_sf(coords = c("x", "y"), crs= sf::st_crs(data)) |>
    dplyr::group_by(line_id) |>
    #see Edzer's answer here:https://github.com/r-spatial/sf/issues/851
    #do_union=FALSE is needed.
    dplyr::summarize(do_union = FALSE) |>
    sf::st_cast("LINESTRING") |> # cast linestring type
    sf::st_as_sf() |>
    dplyr::ungroup()

  return(trajectories_out)
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

# testthat::test_that("exposure_LS normalize range groups", {
#   LS_test = exposure_LS(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "range",
#                         time_data = dateTime, bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
#   LS =  test_exposure_LS(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "range",
#                          time_data = dateTime, bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
#   names(LS) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
#                 "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
#   testthat::expect_equal(LS_test,LS)
# })

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

