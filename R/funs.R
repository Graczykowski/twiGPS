
# handling crs and selecting days from data

start_processing = function(data, x, y, NA_val, env_data, grid_extent,
                            input_crs, output_crs){

  # get spatial data
  if (all(class(data) == "data.frame")) {
    if (missing(input_crs)){ # implement so missing works
      input_crs = ""
    }
    if (!missing(NA_val)) {
      message(paste0("Removing rows containing default NA value ", NA_val, " in x and y columns"))
      data = data[data[x] != NA_val & data[y] != NA_val,]
    }
    if (any(is.na(data[,c(x,y)]))){
      message("Removing rows containing NA in x and y columns")
      data = tidyr::drop_na(data, x, y)
    }
    data_points = terra::vect(x = data, geom = c(x, y), crs = input_crs)

  } else if (any(class(data) == "sf")){
    data_points = terra::vect(data)

    if (!missing(input_crs)){
      message("Ignoring input_crs argument")
    }

  } else if (any(class(data) == "SpatVector")) {
    data_points = data
    if (!missing(input_crs)){
      message("Ignoring input_crs argument")
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
      data_proj = terra::crop(data_proj, terra::ext(grid_extent))
    }
  }

  return(data_proj)

}

calc_grid = function(x, bandwidth, cellsize, env_data, grid_extent, is_LS = FALSE){

  if (!missing(grid_extent)){

    extent = grid_extent
  } else {
    extent = terra::ext(x)
    if (!missing(bandwidth)) {
      extent = c(terra::xmin(extent) - bandwidth, terra::xmax(extent) + bandwidth,
                 terra::ymin(extent) - bandwidth, terra::ymax(extent) + bandwidth)
    }
  }


  if(!missing(cellsize) && is.numeric(cellsize) && cellsize > 0) { # cellsize included

    if (suppressWarnings(!is.na(terra::is.lonlat(x)) && terra::is.lonlat(x))){
      # crs units in degrees
      # is.na if empty crs to skip error
      warning("Cellsize is not stable - cells are not rectangular")

      if (!missing(bandwidth) && bandwidth > 0.1 && is_LS == FALSE) {
        message(paste0("CRS is in lontitude/latitude and bandwidth is ", bandwidth,
      ".Bandwidth is calculated in CRS units. Is bandwidth in correct unit?"))
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
    warning("Cellsize not provided - using cellsize from env_data")

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
    stop("Missing cellsize - provide cellsize or env_data")
  }
  return(grid_rast)
}

# transform env_data SpatVector to env_data SpatRaster

env_vect = function(env, env_buff, env_field, grid){

  if (!missing(env_buff)) {
    if (is.numeric(env_buff) & env_buff > 0){
      env = terra::buffer(env, env_buff)
    } else {
      warning("Incorrect env_buff argument. Creation of buffer skipped.")
    }
  }

  if (!missing(env_field)){
    env_f_enq = rlang::enquo(env_field)
    env = terra::rasterize(env, grid, field = rlang::quo_name(env_f_enq), fun = "sum")
  } else {
    env = terra::rasterize(env, grid, fun = "sum")
  }

  return(env)
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

spat_kde = function(x, ref, bw){

  # coords of each ref cell center
  gx = terra::crds(ref)[,1] |> unique()
  gy = terra::crds(ref)[,2] |> unique()

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
  ax <- outer(ref_uq_x, terra::crds(points)[, 1], "-" ) ^ 2
  ay <- outer(ref_uq_y, terra::crds(points)[, 2], "-" ) ^ 2

  # points within quadratic search radius
  ax_T = ax <= bw ^ 2
  ay_T = ay <= bw ^ 2

  # every x row index
  positions = matrix(1:nrow(ax), ncol = 1)

  # calculate KDE for each column seperately

  density_mx = apply(positions,  1, function(xrow) {

    # which point's x coord within possible search radius
    cols_T = which(ax_T[xrow,] == TRUE)

    # matrix cell coordinates of points within quadratic search radius
    rows_T = which(ay_T[,cols_T] == TRUE, arr.ind = TRUE)

    if (length(cols_T) <= 1){ # if only one or zero points within quadratic search radius
      suppressWarnings({rows_T = rows_T |> array(dim = c(length(rows_T), 1)) |> cbind(1)})
      colnames(rows_T) = c("row", "col")
    }
    # calculate distance within quadratic search radius
    df_row = rows_T |> as.data.frame() |>
      dplyr::mutate(sum_val = ay[cbind(row, cols_T[col])] + ax[cbind(xrow, cols_T[col])]) |>
      dplyr::filter(sum_val <= bw ^ 2) |> # filter points within search radius based on real distance
      dplyr::mutate(sum_val = (sum_val / bw ^ 2 * (-1) + 1) ^ 2) # calculate impact of each point on cells

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
    #length betweeen each point
    ax_p = (vect_x - x) ^ 2
    ay_p = (vect_y - y) ^ 2

    # boolean if dist within search radius
    ax_p_T = ax_p <= bw ^ 2
    ay_p_T = ay_p <= bw ^ 2

    # choose dist^ 2 within search radius
    x_T = which(ax_p_T == TRUE)
    y_T = which(ay_p_T == TRUE)

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
  gx = terra::crds(ref)[,1] |> unique()
  gy = terra::crds(ref)[,2] |> unique()

  # distance from each x/y ref cell center to each x/y x points (squared)

  kde_val = kde(x, gx, gy, bw)
  dr_val = kde_points(terra::crds(x)[,1], terra::crds(x)[,2], bw)
  # DR EVAL
  dr_calc = ecdf(dr_val)(as.vector(kde_val))

  output_list = list(x = gx, y = gy, z = dr_calc)

  # create df of x/y coords and val
  pts = data.frame(expand.grid(x = output_list$x, y = output_list$y),
                   z = output_list$z)
  # create points SpatVector
  pts = terra::vect(pts, geom = c("x", "y"))
  # rasterize to ref
  return(terra::rasterize(pts, ref, field = "z"))
}
