source("testthat-helper.R")
testthat::test_that("exposure_KDE normalize range", {
  KDE_test =  exposure_KDE(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "range",
                        bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  suppressWarnings({KDE =  test_exposure_KDE(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "range",
                         bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")})
  names(KDE) = "activity_space"
  testthat::expect_true(terra::all.equal(KDE_test,KDE))
  })

testthat::test_that("exposure_KDE normalize center", {
  KDE_test =  exposure_KDE(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "center",
                                bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  suppressWarnings({KDE =  test_exposure_KDE(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "center",
                         bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")})
  names(KDE) = "activity_space"
  testthat::expect_true(terra::all.equal(KDE_test,KDE))
})

testthat::test_that("exposure_KDE normalize scale", {
  KDE_test =  exposure_KDE(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "scale",
                                bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  suppressWarnings({KDE =  test_exposure_KDE(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "scale",
                         bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")})
  names(KDE) = "activity_space"
  testthat::expect_true(terra::all.equal(KDE_test,KDE))
})

testthat::test_that("exposure_KDE normalize standardize", {
  KDE_test =  exposure_KDE(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "standardize",
                                bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  suppressWarnings({KDE =  test_exposure_KDE(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "standardize",
                         bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")})
  names(KDE) = "activity_space"
  testthat::expect_true(terra::all.equal(KDE_test,KDE))
})

testthat::test_that("exposure_KDE normalize range groups", {
  KDE_test = exposure_KDE(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "range",
                        bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  suppressWarnings({KDE =  test_exposure_KDE(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "range",
                         bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")})
  names(KDE) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
                "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
  testthat::expect_true(terra::all.equal(KDE_test,KDE))
})

testthat::test_that("exposure_KDE normalize range groups", {
  KDE_test = exposure_KDE(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "range",
                        bandwidth = 200, norm_group = TRUE, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  suppressWarnings({KDE =  test_exposure_KDE(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "range",
                         bandwidth = 200, norm_group = TRUE, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")})
  names(KDE) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
                 "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
  testthat::expect_true(terra::all.equal(KDE_test,KDE))
})

testthat::test_that("exposure_KDE normalize center groups", {
  KDE_test =  exposure_KDE(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "center",
                                bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  suppressWarnings({KDE =  test_exposure_KDE(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "center",
                         bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")})
  names(KDE) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
                 "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
  testthat::expect_true(terra::all.equal(KDE_test,KDE))
})

testthat::test_that("exposure_KDE normalize scale groups", {
  KDE_test =  exposure_KDE(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "scale",
                                bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  suppressWarnings({KDE =  test_exposure_KDE(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "scale",
                         bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")})
  names(KDE) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
                 "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
  testthat::expect_true(terra::all.equal(KDE_test,KDE))
})

testthat::test_that("exposure_KDE normalize standardize groups", {
  KDE_test =  exposure_KDE(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "standardize",
                                bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  suppressWarnings({KDE =  test_exposure_KDE(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "standardize",
                         bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")})
  names(KDE) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
                 "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
  testthat::expect_true(terra::all.equal(KDE_test,KDE))
})

testthat::test_that("exposure_KDE env_data", {
  ndvi_data = terra::rast(system.file("extdata/landsat_ndvi.tif", package = "twiGPS"))

  KDE_test =  exposure_KDE(data = geolife_sandiego ,coords = c("lon", "lat"), cellsize = 50,
                                bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611", env_data = ndvi_data)
  suppressWarnings({KDE =  test_exposure_KDE(data = geolife_sandiego, x = lon, y = lat, cellsize = 50,
                         bandwidth = 200, env_data = ndvi_data, input_crs = "EPSG:4326", output_crs = "EPSG:32611")})
  names(KDE) = "env_exposure"

  testthat::expect_true(terra::all.equal(KDE_test,KDE))
})
