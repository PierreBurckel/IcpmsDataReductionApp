source("C:/Users/pierr/Desktop/IPGP/R/ICP-MS_process/R/functions.R")

# Setup -------------------------------------------------------------------

element_names <- c("Li", "Be", "B")

matrix1_count_value <- matrix(c(4, 6, 5, 3, 0.75, 1.5, 4, 8, 2), nrow = 3, ncol = 3)
matrix1_count_sd <- matrix1_count_value
matrix1_concentration_value <- matrix(c(8, 12, 10, 4.5, 1.125, 2.25, 1, 2, 0.5), nrow = 3, ncol = 3)
matrix1_concentration_sd <- matrix1_concentration_value
background_equivalent_concentration <- matrix1_concentration_value
background_equivalent_sd <- matrix1_concentration_sd

eudc_count <- EstimationUncertaintyDataCouple$new(element_names, matrix1_count_value, matrix1_count_sd, "sd")
eudc_concentration <- EstimationUncertaintyDataCouple$new(element_names, matrix1_concentration_value, matrix1_concentration_sd, "sd")
eudc_blank <- EstimationUncertaintyDataCouple$new(element_names, matrix1_count_value, matrix1_count_sd, "sd")
eudc_detection_limit <- EstimationUncertaintyDataCouple$new(element_names, background_equivalent_sd * 3.3, matrix(0, nrow = 3, ncol = 3), "sd")
eudc_quantification_limit <- EstimationUncertaintyDataCouple$new(element_names, background_equivalent_sd * 10, matrix(0, nrow = 3, ncol = 3), "sd")
zero_eudc <- EstimationUncertaintyDataCouple$new(element_names, matrix(0, nrow = 3, ncol = 3), matrix(0, nrow = 3, ncol = 3), "sd")

calibration_coefficients <- c(2, 1.5, 0.25)

# Testing -----------------------------------------------------------------

test_that("calculateConcentration function works", {
  expect_equal(calculateConcentration(eudc_count, calibration_coefficients),
               eudc_concentration
                   )
  expect_error(calculateConcentration(matrix1_count_value, calibration_coefficients)
  )
})

test_that("calculateLimit function works", {
  expect_equal(calculateLimit(eudc_blank, "Detection Limit", calibration_coefficients),
               eudc_detection_limit
  )
  expect_equal(calculateLimit(eudc_blank, "Quantification Limit", calibration_coefficients),
               eudc_quantification_limit
  )
  # error because "Quantific Limit" isn't allowed
  expect_error(calculateLimit(eudc_blank, "Quantific Limit", calibration_coefficients)
  )
  # type error, matrix not allowed as first argument
  expect_error(calculateLimit(matrix1_count_value, "Quantification Limit", calibration_coefficients)
  )
  # error because size of columns from eudc_blank and calibration coefficients isn't the same
  expect_error(calculateLimit(eudc_blank, "Quantification Limit", 0)
  )
})

test_that("createBlankEudcFromTemplate function works", {
  expect_equal(createBlankEudcFromTemplate(eudc_count),
               zero_eudc
  )
  expect_error(createBlankEudcFromTemplate(matrix1_count_value)
  )
})
