# Setup -------------------------------------------------------------------

nametotest1 <- "7  Li  [ No Gas ] "
nametotest2 <- "75  As  [ He ] "
nametotest3 <- "238  U  [ No Gas ] "
nametotest4 <- c(nametotest1, nametotest2, nametotest3)
expected1 <- "Li"
expected2 <- "As"
expected3 <- "U"
expected4 <- c(expected1, expected2, expected3)

# Testing -----------------------------------------------------------------

test_that("getElementName function works", {
  expect_identical(getElementName(nametotest1, "Agilent7900"),
                   expected1
                   )
  expect_identical(getElementName(nametotest2, "Agilent7900"),
                   expected2
  )
  expect_identical(getElementName(nametotest3, "Agilent7900"),
                   expected3
  )
  expect_identical(getElementName(nametotest4, "Agilent7900"),
                   expected4
  )
  expect_error(getElementName(nametotest1, "Agilent79000"))
})
