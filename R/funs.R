
# handling crs and selecting days from data

start_processing = function(x, day=NULL, env_data = NULL, data_extent = NULL,
                            start_crs = "WGS84", end_crs = NULL){

  # get spatial data
  if(is.null(day)){
    x_points = terra::vect(x = x, geom = c("lon", "lat"), crs = start_crs)
  }
  else{
    x_points = x |> dplyr::filter(wearDate == day) |>
      terra::vect(geom = c("lon", "lat"), crs = start_crs)

  }

  #clip to extent TODO


  # crs
  if (!is.null(end_crs)){
    x_proj = terra::project(x_points, terra::crs(end_crs))
  } else if (!is.null(env_data)) {
    x_proj = terra::project(x_points, terra::crs(env_data))
  } else {
    x_proj = x_points
  }

  return(x_proj)

}

#raster statistics

rast_stats = function(raster, stats){
  vals = c()
  for (statistic in stats){
    # TODO range outputs 2 values and it moves names of stats output
    if (statistic == "range"){
      range = terra::global(raster, fun = statistic, na.rm = TRUE)
      vals = append(vals, range[,2] - range[,1])
    } else if (statistic == "count"){
      count = raster |> terra::freq() |> dplyr::summarise(n = sum(count)) |>
        as.integer()
      vals = append(vals, freq)
    } else if (statistic == "area") {
      area = raster |> terra::expanse() |> dplyr::select(area) |> as.integer()
      vals = append(vals, area)
    } else {
      vals = append(vals, terra::global(raster, fun = statistic, na.rm = TRUE))
    }


  }

  data = as.data.frame(t(vals))
  colnames(data) = stats
  return(data)

}


# Output raster and statistics calculation

output_calc = function(act_rast, env_rast = NULL, stats = NULL,
                       act_and_env = FALSE){ # TODO act_and_env

  output = list()

  if (is.null(env_rast)){

    output$rast = act_rast

    if (!is.null(stats)){ # activity stats

      activity_stats = rast_stats(act_rast, stats)
      output$stats = activity_stats
    }
  } else{

    output$rast = env_rast

    if (!is.null(stats)){ # env stats

      env_stats = rast_stats(env_rast, stats)
      output$stats = env_stats
    }
  }

  return(output)

}

#Line segment trajectories

trajectories_fun = function(data){
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
      names_to = c(".value", "time"),
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

# simple DR for projected DR



DR_simple = function(data, h, kernel="Gaussian",
                     xlim=NULL, ylim=NULL, CL_lev = NULL, x_res=201, y_res = 201, ...){
  if(is.null(h)){
    cat("Please select smoothing bandwidth.")
    break
  }
  if(ncol(data)>2){
    cat("The current version only support 2D.")
    break
  }
  if(is.null(xlim)){
    xlim = range(data[,1])
  }
  if(is.null(ylim)){
    ylim = range(data[,2])
  }
  n = nrow(data)

  ### derived parameters
  x_seq = seq(from=xlim[1], to=xlim[2], length.out=x_res)
  y_seq = seq(from=ylim[1], to=ylim[2], length.out=y_res)
  gr0 = expand.grid(x_seq,y_seq)

  ### alpha value for each point
  if(kernel=="Gaussian"){
    D_kde = TDA::kde(X = data,Grid = data,h = h)
    gr0_kde = TDA::kde(X = data,Grid = gr0,h = h)
  }
  if(kernel=="uniform"){
    D_kde = rowSums(D_nn$nn.idx!=0)/n
    gr0_nn = RANN::nn2(data,query = gr0,k = n,searchtype = "radius",radius = h)
    gr0_kde = rowSums(gr0_nn$nn.idx!=0)/n
  }
  if(kernel=="quartic"){
    D_nn = RANN::nn2(data,data,k = n,searchtype = "radius",radius = h)
    D_val = (1-D_nn$nn.dists^2/h^2)^2*(D_nn$nn.idx!=0)
    D_kde = rowSums(D_val,na.rm = T)
    gr0_nn = RANN::nn2(data,query = gr0,k = n,searchtype = "radius",radius = h)
    gr0_val = (1-gr0_nn$nn.dists^2/h^2)^2*(gr0_nn$nn.idx!=0)
    gr0_kde = rowSums(gr0_val,na.rm = T)
  }
  D_alpha = rank(D_kde)/n
  gr0_alpha = ecdf(D_kde)(gr0_kde)


  out_put = list()
  out_put$h = h
  out_put$x_grid = x_seq
  out_put$y_grid = y_seq
  out_put$gr_alpha = gr0_alpha
  out_put$data_alpha = D_alpha
  return(out_put)
}
