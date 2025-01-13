
testthat::test_that("exposure_DR normalize range", {
 DR_test =  exposure_DR(data = geolife_sandiego ,coords = c("lon", "lat"), cellsize = 50, normalize = "range",
                        bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
 DR =  test_exposure_DR(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "range",
                        bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
 names(DR) = "activity_space"
  testthat::expect_true(terra::all.equal(DR_test,DR))
})

testthat::test_that("exposure_DR normalize center", {
 DR_test = exposure_DR(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "center",
                               bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
 DR =  test_exposure_DR(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "center",
                        bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
 names(DR) = "activity_space"
  testthat::expect_true(terra::all.equal(DR_test,DR))
})

testthat::test_that("exposure_DR normalize scale", {
 DR_test = exposure_DR(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "scale",
                               bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
 DR =  test_exposure_DR(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "scale",
                        bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
 names(DR) = "activity_space"
  testthat::expect_true(terra::all.equal(DR_test,DR))
})

testthat::test_that("exposure_DR normalize standardize", {
 DR_test = exposure_DR(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "standardize",
                               bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
 DR =  test_exposure_DR(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "standardize",
                        bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
 names(DR) = "activity_space"

  testthat::expect_true(terra::all.equal(DR_test,DR))
})


testthat::test_that("exposure_DR normalize range groups", {
 DR_test = exposure_DR(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "range",
                       bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
 DR =  test_exposure_DR(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "range",
                        bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
 names(DR) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
               "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
 testthat::expect_true(terra::all.equal(DR_test,DR))
})

testthat::test_that("exposure_DR normalize range groups", {
 DR_test = exposure_DR(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "range",
                       bandwidth = 200, norm_group = TRUE, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
 DR =  test_exposure_DR(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "range",
                        bandwidth = 200, norm_group = TRUE, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
 names(DR) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
               "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
 testthat::expect_true(terra::all.equal(DR_test,DR))
})

testthat::test_that("exposure_DR normalize center groups", {
 DR_test = exposure_DR(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "center",
                               bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
 DR =  test_exposure_DR(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "center",
                        bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
 names(DR) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
               "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
 testthat::expect_true(terra::all.equal(DR_test,DR))
})

testthat::test_that("exposure_DR normalize scale groups", {
 DR_test = exposure_DR(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "scale",
                               bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
 DR =  test_exposure_DR(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "scale",
                        bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
 names(DR) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
               "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
   testthat::expect_true(terra::all.equal(DR_test,DR))
})

testthat::test_that("exposure_DR normalize standardize groups", {
 DR_test = exposure_DR(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50, normalize = "standardize",
                               bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
 DR =  test_exposure_DR(data = geolife_sandiego, x = lon, y = lat, cellsize = 50, normalize = TRUE, norm_method = "standardize",
                        bandwidth = 200, group_split = date, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
 names(DR) = c("2011-08-17", "2011-08-18", "2011-08-19", "2011-08-20",
               "2011-08-21", "2011-08-22", "2011-08-23", "2011-08-24")
  testthat::expect_true(terra::all.equal(DR_test,DR))
})

testthat::test_that("exposure_DR env_data", {
  ndvi_data = terra::rast(system.file("extdata/landsat_ndvi.tif", package = "twiGPS"))

 DR_test = exposure_DR(data = geolife_sandiego, coords = c("lon", "lat"), cellsize = 50,
                               bandwidth = 200, input_crs = "EPSG:4326", output_crs = "EPSG:32611", env_data = ndvi_data)
 DR =  test_exposure_DR(data = geolife_sandiego, x = lon, y = lat, cellsize = 50,
                        bandwidth = 200, env_data = ndvi_data, input_crs = "EPSG:4326", output_crs = "EPSG:32611")
 names(DR) = "env_exposure"
  testthat::expect_true(terra::all.equal(DR_test,DR))
})
