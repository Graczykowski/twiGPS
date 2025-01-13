
testthat::test_that("exposure_PO normalize range", {
  PO_test = exposure_PO(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "range",
                   input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  PO =  test_exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "range",
                    input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(PO) = "activity_space"
  testthat::expect_true(terra::all.equal(PO_test,PO))
})

testthat::test_that("exposure_PO normalize center", {
  PO_test = exposure_PO(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "center",
                        input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  PO =  test_exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "center",
                         input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(PO) = "activity_space"
  testthat::expect_true(terra::all.equal(PO_test,PO))
})

testthat::test_that("exposure_PO normalize scale", {
  PO_test = exposure_PO(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "scale",
                        input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  PO =  test_exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "scale",
                         input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(PO) = "activity_space"
  testthat::expect_true(terra::all.equal(PO_test,PO))
})

testthat::test_that("exposure_PO normalize standardize", {
  PO_test = exposure_PO(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "standardize",
                        input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  PO =  test_exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "standardize",
                         input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(PO) = "activity_space"
  testthat::expect_true(terra::all.equal(PO_test,PO))
})


testthat::test_that("exposure_PO normalize range groups", {
  PO_test = exposure_PO(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "range",
                        group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  PO =  test_exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "range",
                         group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(PO) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
                "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
  testthat::expect_true(terra::all.equal(PO_test,PO))
})

testthat::test_that("exposure_PO normalize range groups", {
  PO_test = exposure_PO(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "range",
                        norm_group = TRUE, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  PO =  test_exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "range",
                         norm_group = TRUE, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(PO) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
                "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
  testthat::expect_true(terra::all.equal(PO_test,PO))
})

testthat::test_that("exposure_PO normalize center groups", {
  PO_test = exposure_PO(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "center",
                                group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  PO =  test_exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "center",
                         group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(PO) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
                "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
  testthat::expect_true(terra::all.equal(PO_test,PO))
})

testthat::test_that("exposure_PO normalize scale groups", {
  PO_test = exposure_PO(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "scale",
                                group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  PO =  test_exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "scale",
                         group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(PO) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
                "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
  testthat::expect_true(terra::all.equal(PO_test,PO))
})

testthat::test_that("exposure_PO normalize standardize groups", {
  PO_test = exposure_PO(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "standardize",
                                group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  PO =  test_exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "standardize",
                         group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(PO) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
                "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
  testthat::expect_true(terra::all.equal(PO_test,PO))
})

testthat::test_that("exposure_PO env_data", {
  ndvi_data = terra::rast(system.file("extdata/landsat_ndvi.tif", package = "twiGPS"))

  PO_test = twiGPS::exposure_PO(data = twiGPS::geolife_sandiego, coords = c("lon", "lat"), cellsize = 50,
                                 input_crs = "EPSG:4326", output_crs = "EPSG:32611", env_data = ndvi_data)
  PO =  test_exposure_PO(data = geolife_sandiego, x = lon, y = lat, cellsize = 50,
                         env_data = ndvi_data, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
  names(PO) = "env_exposure"
  testthat::expect_true(terra::all.equal(PO_test,PO))
})
