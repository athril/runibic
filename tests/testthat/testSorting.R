context("Sorting using indexes")

test_that("Sorting based on indexes: unisort", {
  A <- matrix(c(15,3,10,11,12,10,1,7,10,2,6,11,8,6,16,9), nrow = 2, byrow = TRUE)
  result <- matrix(c(6,1,7,2,5,3,4,0,1,2,5,4,7,0,3,6), nrow = 2, byrow = TRUE)

  B <- unisort(A)
  expect_that( B, is_a("matrix"))
  expect_that( B, equals(result))
})

