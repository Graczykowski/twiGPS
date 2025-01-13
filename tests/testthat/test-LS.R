# library(twiGPS)

testthat::test_that("exposure_LS normalize range", {
  LS_test = exposure_LS(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "range",
                        time_data = dateTime, bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")

  LS =  test_exposure_LS(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "range",
                         time_data = dateTime, bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(LS) = "activity_space"
  testthat::expect_true(terra::all.equal(LS_test,LS))
})

testthat::test_that("exposure_LS normalize center", {
  LS_test =  exposure_LS(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "center",
                                time_data = dateTime, bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  LS =  test_exposure_LS(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "center",
                         time_data = dateTime, bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(LS) = "activity_space"
  testthat::expect_true(terra::all.equal(LS_test,LS))
})

testthat::test_that("exposure_LS normalize scale", {
  LS_test =  exposure_LS(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "scale",
                                time_data = dateTime, bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  LS =  test_exposure_LS(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "scale",
                         time_data = dateTime, bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(LS) = "activity_space"
  testthat::expect_true(terra::all.equal(LS_test,LS))
})

testthat::test_that("exposure_LS normalize standardize", {
  LS_test =  exposure_LS(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "standardize",
                                time_data = dateTime, bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  LS =  test_exposure_LS(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "standardize",
                         time_data = dateTime, bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(LS) = "activity_space"
  testthat::expect_true(terra::all.equal(LS_test,LS))
})

testthat::test_that("exposure_LS normalize range groups", {
  LS_test = exposure_LS(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "range",
                        time_data = dateTime, bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  LS =  test_exposure_LS(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "range",
                         time_data = dateTime, bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(LS) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
                "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
  testthat::expect_true(terra::all.equal(LS_test,LS))
})

testthat::test_that("exposure_LS normalize range norm groups", {
  LS_test = exposure_LS(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "range",
                        time_data = dateTime, bandwidth = 200, norm_group = TRUE, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  LS =  test_exposure_LS(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "range",
                         time_data = dateTime, bandwidth = 200, norm_group = TRUE, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(LS) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
                "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
  testthat::expect_true(terra::all.equal(LS_test,LS))
})

testthat::test_that("exposure_LS normalize center groups", {
  LS_test =  exposure_LS(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "center",
                                time_data = dateTime, bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  LS =  test_exposure_LS(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "center",
                         time_data = dateTime, bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(LS) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
                "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
  testthat::expect_true(terra::all.equal(LS_test,LS))
})

testthat::test_that("exposure_LS normalize scale groups", {
  LS_test =  exposure_LS(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "scale",
                                time_data = dateTime, bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  LS =  test_exposure_LS(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "scale",
                         time_data = dateTime, bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(LS) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
                "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
  testthat::expect_true(terra::all.equal(LS_test,LS))
})

testthat::test_that("exposure_LS normalize standardize groups", {
  LS_test =  exposure_LS(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "standardize",
                                time_data = dateTime, bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  LS =  test_exposure_LS(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "standardize",
                         time_data = dateTime, bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(LS) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
                "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
  testthat::expect_true(terra::all.equal(LS_test,LS))
})

testthat::test_that("exposure_LS env_data", {
  ndvi_data = terra::rast(system.file("extdata/landsat_ndvi.tif", package = "twiGPS"))

  LS_test =  exposure_LS(data = geolife_sandiego ,coords = c("lon", "lat"), cellsize = 50,
                                time_data = dateTime, bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611", env_data = ndvi_data)
  LS =  test_exposure_LS(data = geolife_sandiego, x = lon, y = lat, cellsize = 50,
                         time_data = dateTime, bandwidth = 200, env_data = ndvi_data, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(LS) = "env_exposure"
  testthat::expect_true(terra::all.equal(LS_test,LS))
})

