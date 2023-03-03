# Setup -------------------------------------------------------------------

element_names <- c("Li", "Be", "B")

matrix1_value <- matrix(seq(9), nrow = 3, ncol = 3)
matrix1_sd <- matrix(rep(seq(3), 3), nrow = 3, ncol = 3)
matrix2_value <- matrix(11:19, nrow = 3, ncol = 3)
matrix2_sd <- matrix(rep(-1, 9), nrow = 3, ncol = 3)
matrix_zero <- matrix(0, nrow = 3, ncol = 3)
matrix1_rsd <- matrix1_sd / matrix1_value * 100
matrix3_value <- matrix1_value
matrix3_sd <- matrix1_sd
matrix3_value[,2] <- matrix2_value[,2]
matrix3_sd[,2] <- matrix2_sd[,2]

colnames(matrix1_value) <- element_names
colnames(matrix1_sd) <- element_names
colnames(matrix2_value) <- element_names
colnames(matrix2_sd) <- element_names
colnames(matrix_zero) <- element_names
colnames(matrix1_rsd) <- element_names

matrix_diff_value <- matrix1_value - matrix2_value
matrix_sum_value <- matrix1_value + matrix2_value
matrix_mult_value <- matrix1_value * matrix2_value
matrix_div_value <- matrix1_value / matrix2_value

matrix_diff_sd <- sqrt(matrix1_sd^2 + matrix2_sd^2)
matrix_sum_sd <- sqrt(matrix1_sd^2 + matrix2_sd^2)
matrix_mult_sd <- sqrt((matrix1_sd/matrix1_value)^2 + (matrix2_sd/matrix2_value)^2)*matrix_mult_value
matrix_div_sd <- sqrt((matrix1_sd/matrix1_value)^2 + (matrix2_sd/matrix2_value)^2)*matrix_div_value

eudc1 <- EstimationUncertaintyDataCouple$new(element_names, matrix1_value, matrix1_sd, "sd")
eudc2 <- EstimationUncertaintyDataCouple$new(element_names, matrix2_value, matrix2_sd, "sd")
eudc3 <- eudc1
eudc3[ ,2] <- eudc2[ ,2]

eudc_diff <- eudc1$subtract(eudc2)
eudc_sum <- eudc1$add(eudc2)
eudc_mult <- eudc1$multiplyBy(eudc2)
eudc_div <- eudc1$divideBy(eudc2)

# Testing -----------------------------------------------------------------

test_that("EUDC methods work", {
  expect_equal(list(matrix_diff_value, matrix_diff_sd),
               list(eudc_diff$getEstimation(), eudc_diff$getSd())
                   )
  expect_equal(list(matrix_sum_value, matrix_sum_sd),
               list(eudc_sum$getEstimation(), eudc_sum$getSd())
  )
  # expect_equal is more tolerant on rounding
  expect_equal(list(matrix_mult_value, matrix_mult_sd),
                   list(eudc_mult$getEstimation(), eudc_mult$getSd())
  )
  # expect_equal is more tolerant on rounding
  expect_equal(list(matrix_div_value, matrix_div_sd),
                   list(eudc_div$getEstimation(), eudc_div$getSd())
  )
  expect_equal(list(matrix_zero, matrix_zero),
               list(eudc_diff$zeroOutNegatives()$getEstimation(), eudc_diff$zeroOutNegatives()$getSd())
  )
  expect_identical(eudc_diff$getElementFullNames(),
                   eudc_sum$getElementFullNames(),
                   eudc_mult$getElementFullNames(),
                   eudc_div$getElementFullNames(),
                   element_names
  )
  expect_identical(eudc1$getElementFullNames(),
                   eudc2$getElementFullNames(),
                   element_names
  )
  expect_equal(eudc1$getSampleNumber(),
               3
  )
  expect_equal(eudc1$getElementNumber(),
               3
  )
  expect_equal(eudc1$getRsd(),
               matrix1_rsd
  )
  expect_identical(eudc1$createEmptyMatrixOfEudcSize(),
                   matrix(nrow = 3, ncol = 3)
  )
  expect_equal(eudc3,
                   EstimationUncertaintyDataCouple$new(element_names, matrix3_value, matrix3_sd, "sd")
  )
  expect_error(EstimationUncertaintyDataCouple$new(element_names, matrix1_value, matrix1_sd, "rrsd"))
})
