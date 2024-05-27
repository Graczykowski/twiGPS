#' Geolife GPS Trajectory 1.3 Data San Diego Subset
#'
#' @description
#' Subset of GPS real-life time-stamped dataset from 1 participant collected in Geolife GPS Trajectory 1.3 project
#'  between 17th August 2011 and 20th August 2011. Area of subset is limited to center of San Diego city.
#' @format
#' A data frame with 1068 rows and 5 columns:
#' \describe{
#'    \item{lon, lat}{longtitude and latitude coordinates of GPS points}
#'    \item{date}{date when GPS point was measured}
#'    \item{time}{time when GPS point was measured}
#'    \item{dateTime}{time-stamp when GPS point was measured}
#' }
#' @references Yu Zheng, Wei-Ying Ma, Xing Xie. 2010. “GeoLife: A Collaborative Social Networking Service Among User, Location and Trajectory.” IEEE Data Engineering Bulletin 33 (2): 32–44. https://www.microsoft.com/en-us/research/wp-content/uploads/2016/02/GeoLife-A20Collaborative20Social20Networking20Service20among20User2C20Location20and20Trajectory.pdf.
#'
#' Zheng, Yu, Quannan Li, Yukun Chen, Xing Xie, and Wei-Ying Ma. 2008. “Understanding Mobility Based on GPS Data.” Proceedings of the 10th International Conference on Ubiquitous Computing, September. https://doi.org/10.1145/1409635.1409677.
#'
#' Zheng, Yu, Lizhu Zhang, Xing Xie, and Wei-Ying Ma. 2009. “Mining Interesting Locations and Travel Sequences from GPS Trajectories.” Proceedings of the 18th International Conference on World Wide Web, April. https://doi.org/10.1145/1526709.1526816.
"geolife_sandiego"

#' Landsat 9 NDVI image San Diego subset
#'
#' @description
#' San Diego subset of raster NDVI image calculated from Landsat 9 scene
#' LC90400372024053LGN00 acquired on 22nd February 2024. Raster is in "EPSG:32611" or "WGS84/UTM zone 11N" CRS.
#'
#' @usage data(landsat_ndvi)
#' @format RasterLayer with 390 rows, 304 columns and 1 layer
#'
#' @references U.S. Geological Survey
"landsat_ndvi"
