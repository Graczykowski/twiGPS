#' @importFrom rlang .data
# handling crs and selecting days from data

start_processing = function(data, coords, NA_val, env_data, grid_extent,
                            input_crs, output_crs, verbose){

  # get spatial data
  if (inherits(data, "data.frame") && !inherits(data, "sf")) {
    if (missing(input_crs)){ # implement so missing works
      input_crs = ""
    }
    if (!missing(NA_val)) {
      if (is.numeric(NA_val)){
        n_row = nrow(data)
        data = dplyr::filter(data, .data[[coords]] != NA_val)
        data[data[coords[1]] != NA_val & data[coords[2]] != NA_val,] # to improve
        if (verbose) {
          number_row = n_row - nrow(data)
          warning(paste0("Removed ", number_row, ifelse(number_row == 1, " row", " rows"), " containing default NA value ",
                         NA_val, " in 'coords' columns names"))
        }
      } else if (verbose) {
        warning("Argument 'NA_val' should be numeric. Ignoring 'NA_val' argument")
      }

    }
    if (any(is.na(data[,coords]))){
      n_row_NA = nrow(data)
      data = tidyr::drop_na(data, coords)
      if (verbose) {
        n_na_rows = n_row_NA - nrow(data)
        warning(paste0("Removed ", n_na_rows, ifelse(n_na_rows == 1, " row", " rows"),
                       " containing NA in 'coords' columns names"))
      }

    }
    data_points = terra::vect(x = data, geom = coords, crs = input_crs)

  } else if (inherits(data, "sf") && all(sf::st_geometry_type(data) == "POINT")){
    data_points = terra::vect(data)

    if (!missing(input_crs) && verbose){
      warning("Argument 'input_crs' is ignored")
    }

  } else if (inherits(data, "SpatVector") && terra::geomtype(data) == "points") {
    data_points = data
    if (!missing(input_crs) && verbose){
      warning("Argument 'input_crs' is ignored")
    }
  } else {
    stop("Argument 'data' should be data.frame, sf or SpatVector class")
  }

  data_crs = terra::crs(data_points)

  # crs
  if (!missing(output_crs)){
    if (data_crs == "") {
      if (verbose) {
        message("Setting 'data' CRS to 'output_crs'")
      }
      terra::crs(data_points) = output_crs
      data_proj = data_points
    } else {
      if (verbose) {
        message("Projecting 'data' to 'output_crs'")
      }
      data_proj = terra::project(data_points, output_crs)
    }
  } else if (!missing(grid_extent) && inherits(grid_extent, "SpatRaster")) {
    if (data_crs == "") {
      if (verbose) {
        message("Setting 'data' CRS to 'grid_extent' CRS")
      }
      terra::crs(data_points) = terra::crs(grid_extent)
      data_proj = data_points
    } else {
      if (verbose) {
        message("Projecting 'data' to 'grid_extent' CRS")
      }
      data_proj = terra::project(data_points, grid_extent)
    }
  } else if (data_crs != "") { # any invalid/empty crs
    data_proj = data_points
  } else if (!missing(env_data)) {
    if (verbose) {
      message("Setting 'data' CRS to 'env_data' CRS")
    }
    terra::crs(data_points) = terra::crs(env_data)
    data_proj = data_points
  } else {
    if (verbose) {
      message("CRS is not specified")
    }
    data_proj = data_points
  }
  # crop to grid_extent
  if (!missing(grid_extent)){
    if (inherits(grid_extent, "SpatExtent")){
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

calc_grid = function(x, bandwidth, cellsize, env_data, grid_extent, is_LS = FALSE, verbose){

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
      if (verbose) {
        warning("Cellsize is not stable - cells are not rectangular")
      }

      if (!missing(bandwidth) && bandwidth > 0.1 && !is_LS && verbose) {
        message(paste0("CRS is in lontitude/latitude and 'bandwidth' is ", bandwidth,
      " - 'bandwidth' is calculated in CRS units. Is 'bandwidth' in correct unit?"))
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

  } else if  (!missing(env_data) && inherits(env_data, "SpatRaster")){ #if incorrect cellsize and env_data exists
    if (verbose) {
      warning("Argument 'cellsize' should be positive and single numeric value - cellsize from 'env_data' was used")
    }
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
    stop("Argument 'cellsize' should be positive and single numeric value - valid 'cellsize' or 'env_data' argument should be provided")
  }
  return(grid_rast)
}

# transform env_data SpatVector to env_data SpatRaster

env_vect = function(env, env_buff, env_field, grid, verbose){

  if (!missing(env_buff)) { # create buffer around vector data
    if (length(env_buff) == 1 && is.numeric(env_buff) && env_buff > 0){
      env_proj = terra::project(env, terra::crs(grid))
      ## TERRA::BUFFER SOMETIMES CRUSHES WITH BIG LINE VECTOR DATASETS
      #env = terra::buffer(env_proj, env_buff)
      ## TEMPORARY FIX: USE SF::ST_BUFFER
      env = env_proj |> sf::st_as_sf() |> sf::st_buffer(env_buff) |> terra::vect()

    } else if (verbose) {
      warning("Argument 'env_buff' should be positive and single numeric value - creation of buffer was ignored")
    }
  }

  if (!missing(env_field)){ # choose field to rasterize
    env_f_quo = rlang::quo_name(rlang::enquo(env_field))
    if (env_f_quo %in% terra::names(env)){
      env = terra::rasterize(env, grid, field = env_f_quo, fun = "sum")
    } else {
      if (verbose) {
        warning("Column name in 'env_field' argument is not in 'data'. Argument 'env_field' was ignored")
      }
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


#Line segment trajectories

trajectories_fun = function(data){

  t_name = c(names(data))[grepl('time', c(names(data)))]
  lst_name = paste(t_name[which.max(nchar(t_name))], "_1", sep = "")

  trajectories_out = data |>
    sf::st_as_sf() |> # change to sf object to change points to linestring later
    #filter to the study id defined by the argument
    #define start and end points of line
    dplyr::mutate(
      line_id = dplyr::row_number(),#an id for each "line segment"
      x_start= sf::st_coordinates(.data$geometry)[,1],
      y_start= sf::st_coordinates(.data$geometry)[,2],
      x_end = dplyr::lead(.data$x_start),
      y_end = dplyr::lead(.data$y_start)
    ) |>
    dplyr::ungroup() |>
    sf::st_set_geometry(NULL) |>
    #exclude the last observation, which has no "lead", and will be missing.
    dplyr::filter(is.na(.data$x_end)==FALSE) |>
    # Make the data long form so that each point has two observations
    tidyr::pivot_longer(
      # select variables to pivot longer.
      cols = c(.data$x_start, .data$y_start, .data$x_end, .data$y_end),
      #value goes to "x/y", and time goes to "_start/end"
      names_to = c(".value", lst_name),
      names_repair = "unique",
      names_sep = "_"#the separator for the column name
    ) |>
    # create sf object once again
    sf::st_as_sf(coords = c("x", "y"), crs= sf::st_crs(data)) |>
    dplyr::group_by(.data$line_id) |>
    #see Edzer's answer here:https://github.com/r-spatial/sf/issues/851
    #do_union=FALSE is needed.
    dplyr::summarize(do_union = FALSE) |>
    sf::st_cast("LINESTRING") |> # cast linestring type
    sf::st_as_sf() |>
    dplyr::ungroup()

  return(trajectories_out)
}

spat_kde = function(x, ref, bw){

  ref_coords = terra::crds(ref)
  # coords of each ref cell center
  gx = ref_coords[,1] |> unique()
  gy = ref_coords[,2] |> unique()

  # distance from each x/y ref cell center to each x/y x points (squared)

  kde_val = kde(x, gx, gy, bw)

  output_list = list(x = gx, y = gy, z = kde_val)

  # create df of x/y coords and val
  pts = data.frame(expand.grid(x = output_list$x, y = output_list$y),
                   z = as.vector(array(output_list$z,  length(output_list$z))))
  # create points SpatVector
  pts = terra::vect(pts, geom = c("x", "y"))
  # rasterize to ref
  return(terra::rasterize(pts, ref, field = "z"))
}

kde = function(points, ref_uq_x, ref_uq_y, bw) {

  points_coords = terra::crds(points)

  ax <- outer(ref_uq_x, points_coords[, 1], "-" ) ^ 2
  ay <- outer(ref_uq_y, points_coords[, 2], "-" ) ^ 2

  # points within quadratic search radius
  ax_T = ax <= bw ^ 2
  ay_T = ay <= bw ^ 2

  # every x row index
  positions = matrix(1:nrow(ax), ncol = 1)

  # calculate KDE for each column seperately

  density_mx = apply(positions,  1, function(xcol) {

    # which point's x coord within possible search radius
    cols_T = which(ax_T[xcol,])

    # matrix cell coordinates of points within quadratic search radius
    rows_T = which(ay_T[,cols_T], arr.ind = TRUE)

    if (length(cols_T) <= 1){ # if only one or zero points within quadratic search radius
      suppressWarnings({rows_T = rows_T |> array(dim = c(length(rows_T), 1)) |> cbind(1)})
      colnames(rows_T) = c("row", "col")
    }
    # calculate distance within quadratic search radius
    df_row = rows_T |> as.data.frame() |>
      dplyr::mutate(sum_val = ay[cbind(row, cols_T[col])] + ax[cbind(xcol, cols_T[col])]) |>
      dplyr::filter(.data$sum_val <= bw ^ 2) |> # filter points within search radius based on real distance
      dplyr::mutate(sum_val = (.data$sum_val / bw ^ 2 * (-1) + 1) ^ 2) # calculate impact of each point on cells

    #create empty matrix
    empty_mx = matrix(NA, nrow = nrow(ay), ncol = length(cols_T))
    # insert values based on cell coordinates
    empty_mx[cbind(df_row$row, df_row$col)] = df_row$sum_val

    # sum impact of every point to calculate KDE
    mx_sum = rowSums(empty_mx, na.rm = TRUE)

    return(mx_sum)})

  # transpose matrix and final KDE calculations
  out = t(density_mx) * (3/pi) / bw ^ 2
  return(out)
}

kde_points = function(vect_x, vect_y, bw){

  # apply for each cell coords seperately
  kde_p = mapply(function(x, y) {
    #distance betweeen each point
    ax_p = (vect_x - x) ^ 2
    ay_p = (vect_y - y) ^ 2

    # boolean if dist within search radius
    ax_p_T = ax_p <= bw ^ 2
    ay_p_T = ay_p <= bw ^ 2

    # choose dist^ 2 within search radius
    x_T = which(ax_p_T)
    y_T = which(ay_p_T)

    # index of coords of which both x and y within search radius
    xy_T = intersect(x_T, y_T)

    # calculate dist ^ 2 and choose thoso within search radius
    xy_val = ax_p[xy_T] + ay_p[xy_T]
    xy_val = xy_val[xy_val < bw ^ 2]

    #
    xy_calc = (xy_val / bw ^ 2 * (-1) + 1) ^ 2

    xy_out = sum(xy_calc)


  }, x = vect_x, y =  vect_y) * (3/pi) / bw ^ 2

  return(kde_p)
}

spat_dr = function(x, ref, bw) {
  ref_coords = terra::crds(ref)
  # coords of each ref cell center
  gx = ref_coords[,1] |> unique()
  gy = ref_coords[,2] |> unique()

  # distance from each x/y ref cell center to each x/y x points (squared)

  x_coords = terra::crds(x)

  kde_val = kde(x, gx, gy, bw)
  dr_val = kde_points(x_coords[,1], x_coords[,2], bw)
  # DR EVAL
  dr_calc = stats::ecdf(dr_val)(as.vector(kde_val))

  output_list = list(x = gx, y = gy, z = dr_calc)

  # create df of x/y coords and val
  pts = data.frame(expand.grid(x = output_list$x, y = output_list$y),
                   z = output_list$z)
  # create points SpatVector
  pts = terra::vect(pts, geom = c("x", "y"))
  # rasterize to ref
  return(terra::rasterize(pts, ref, field = "z"))
}
