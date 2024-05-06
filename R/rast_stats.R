#' Raster statistics
#'
#' @description
#' Calculate statistics dataframe from multiple SpatRasters
#'
#' @param ... SpatRaster objects to calculate statistics on
#' @param stats vector of statistics to be calculated. See terra::global() for more details. Added "range", "count" and "area".
#' @param row_names vector the same length as number of ... arguments specifying row names of output data.frame. If NULL row names are ... arguments, values without argument specified will be numbered row-wise.
#'
#' @return Data.frame of length(stats) columns and ... rows
#'
#' @examples
#' statistics = c("count", "area", "min", "max", "range", "mean", "std", 'sum')
#'
#' ndvi_data = terra::rast(system.file("extdata/landsat_ndvi.tif", package = "twsagps"))
#'
#' exposure = exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE,
#'                        input_crs = "EPSG:4326", output_crs = "EPSG:32611")
#' # no names
#' rast_stats(ndvi_data, exposure, stats = statistics)
#'
#' # argument names
#' rast_stats(ndvi = ndvi_data, PO = exposure, stats = statistics)
#'
#' # names
#' rast_stats(ndvi_data, exposure, stats = statistics, row_names = c("NDVI", "PO"))
#'
#' @export

rast_stats =  function(..., stats, row_names = NULL){

  elipsis_list = list(...)
  el_class = sapply(elipsis_list, class)
  if (!all(el_class == "SpatRaster")){
    if (any(el_class == "SpatRaster")){
      warning("Skipping not SpatRaster class arguments")
      elipsis_list = elipsis_list[el_class == "SpatRaster"]
    } else {
      stop("None of arguments is SpatRaster class")
    }
  }

  if (any(sapply(elipsis_list, function(x) {terra::crs(x)}) == "")  && "area" %in% stats) { # when
    warning("Empty crs of one of arguments - area statistic will not be calculated")
    stats = stats[!stats == "area"]
  }
  if (missing(stats) || length(stats) == 0){
    stop("No statistic to be calculated")
  }

  df = data.frame(matrix(ncol = length(stats)))[-1,]
  for (raster in elipsis_list){

    row_df = data.frame(matrix(NA, ncol=1, nrow=1))[-1]
    for (statistic in stats){

      column = switch(
        statistic,
        range = {range = terra::global(raster, fun = "range", na.rm = TRUE)
        val = range[2] - range[1]
        colnames(val) = "range"
        val},
        count = {raster |> terra::freq() |> dplyr::summarise(count = sum(count))},
        area = {raster |> terra::expanse() |> dplyr::select(area)},
        terra::global(raster, fun = statistic, na.rm = TRUE)
      )
      row_df = cbind(row_df, column)
    }
    df = rbind(df, row_df)
  }

  if (is.null(row_names)){
    el_names = names(elipsis_list)
    names_repl = which(el_names == "")
    el_names[names_repl] = names_repl
    rownames(df) = el_names
  } else if (length(row_names) == length(elipsis_list)) {
    rownames(df) = row_names
  } else {
    warning("Incorrect length of row_names argument - default output row names")
  }


  return(df)
}
