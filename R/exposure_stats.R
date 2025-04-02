#' Exposure statistics
#'
#' @description
#' Calculate statistics data.frame from multiple SpatRasters
#'
#'
#' @param ... SpatRaster objects from which statistics are calculated
#' @param stats Character. Single or vector of statistics to be calculated. See [terra::global()] for more details. Added "range", "count" and "area".
# @param layer_names Boolean. Should the layer names be the row names of output data.frame (Default TRUE)
#' @param row_names Character. Single or vector the same length as number of ... arguments specifying row names of output data.frame. If NULL row names are ... arguments, values without argument specified will be numbered row-wise. If ... SpatRasters have more than one layer row_names vector should be a length of sum of layers of all SpatRasters.
# Overrides layer_names argument.
#' @param verbose Boolean. If TRUE amount of output is reduced.
#'
#' @return Data.frame of `length(stats)` columns and `length(...)` rows
#'
#' @examples
#' statistics = c("count", "area", "min", "max", "range", "mean", "std", 'sum')
#'
#' ndvi_data = terra::rast(system.file("extdata/landsat_ndvi.tif",
#'                                      package = "twiGPS"))
#'
#' exposure = exposure_PO(data = geolife_sandiego, coords = c("lon", "lat"),
#'                        cellsize = 50, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
#' # no names
#' exposure_stats(ndvi_data, exposure, stats = statistics)
#'
#' # argument names
#' exposure_stats(ndvi = ndvi_data, PO = exposure, stats = statistics)
#'
#' # names
#' exposure_stats(ndvi_data, exposure, stats = statistics,
#'  row_names = c("NDVI", "Point Exposure"))
#' @seealso [exposure_PO()], [exposure_KDE()], [exposure_DR()], [exposure_LS()]
#' @importFrom rlang .data
#' @export

exposure_stats =  function(..., stats, #layer_names = TRUE,
                           row_names = NULL, verbose = TRUE){

  elipsis_list = list(...)
  el_class = sapply(elipsis_list, class)
  if (!all(el_class == "SpatRaster")){
    if (any(el_class == "SpatRaster")){
      if (verbose){
        warning("Elipsis (...) arguments should be SpatRaster class. Skipping not SpatRaster class arguments")
      }
      elipsis_list = elipsis_list[el_class == "SpatRaster"]
    } else {
      stop("Any of elipsis (...) arguments should be SpatRaster class")
    }
  }

  if (any(sapply(elipsis_list, function(x) {terra::crs(x)}) == "")  && "area" %in% stats) { # when
    if (verbose){
      warning("Empty CRS of one of elipsis (...) arguments - area statistic will not be calculated")
    }
    stats = stats[!stats == "area"]
  }
  if (missing(stats) || length(stats) == 0){
    stop("No statistic to be calculated")
  }

  # if (is.na(as.logical(layer_names))){
  #   stop("Argument 'layer_names' should be interpretable as boolean")
  # }

  rows = sum(sapply(elipsis_list, terra::nlyr))


  df = data.frame(matrix(ncol = length(stats), nrow = rows))
  for (i in 1:length(elipsis_list)){

    r = elipsis_list[[i]]
    r_names = terra::names(r)

    row_df = data.frame(matrix(ncol=length(stats), nrow=terra::nlyr(r)))
    for (n in 1:length(stats)){

      stat = switch(
        stats[n],
        range = {range = terra::global(r, fun = "range", na.rm = TRUE)
          val = range[2] - range[1]
          val},
        count = {terra::global(r, fun = "notNA", na.rm = TRUE)},
        area = {terra::expanse(r)[, "area"]},
        terra::global(r, fun = stats[n], na.rm = TRUE)
      )
      row_df[,n] = stat
    }
    df[i:(terra::nlyr(r) + i-1), ] = row_df
    rownames(df)[i:(terra::nlyr(r) + i-1)] = r_names
  }

  colnames(df) = stats


  if (is.null(row_names)){
    el_names = names(elipsis_list)
    names_repl = which(el_names == "")

    # if (length(elipsis_list) > 1){
    #   el_lay = which(sapply(elipsis_list, terra::nlyr) > 1)
    #
    #   nlay = sapply(elipsis_list[el_lay], FUN =terra::nlyr)
    #
    #   rep_l = sapply(elipsis_list[names_repl], FUN =terra::names)
    #
    # } else {
    #   if (terra::nlyr(elipsis_list[[1]]) > 1){
    #     if(length(names_repl) == 0){
    #       el_names = terra::names(elipsis_list[[1]])
    #     } else {
    #       el_names = paste(el_names, terra::names(elipsis_list[[1]]), sep='@')
    #     }
    #   } else {
    #     el_names[names_repl] = names_repl
    #   }
    # }


    el_names[names_repl] = names_repl
    end_names = c()

    for (l in el_names){
      if (terra::nlyr(elipsis_list[l][[1]]) > 1){
        n_name = terra::names(elipsis_list[l][[1]])
        if (grep('^(?=.)([+-]?([0-9]*)(\\.([0-9]+))?)$', l, perl = TRUE,
                 invert = TRUE, value = TRUE) |> length() > 0) {
          n_name = paste(l, n_name, sep = '@')
        }
        end_names = c(end_names, n_name)
      } else {
        end_names = c(end_names, l)
      }
    }

    rownames(df) = end_names
  } else if (length(row_names) == rows) {
    rownames(df) = row_names
  } else if (verbose){
    if (rows == length(elipsis_list)){
      warning(paste0("Argument 'row_names' (", length(row_names),
                    ") should be same length as number of elipsis (...) arguments (",
                    length(row_names), ") - row names were set to default", sep =''))
    } else {
      warning(paste0("Argument 'row_names' (", length(row_names),
                     ") should be same length as sum of SpatRaster layers of elipsis (...) arguments (",
                     rows, ") - row names were set to default", sep = ''))
    }
  }


  return(df)
}
