context("Calculating LCS")

test_that("Calculating Longest Common Subsequence (LCS): calculateLCS", {
  A = matrix(c(6,1,7,2,5,3,4,0,1,2,5,4,7,0,3,6), nrow=2, byrow=TRUE)
  result=matrix(c(0,0,0,0,0,0,0,0,
           0,0,0,0,0,0,0,0,
           0,1,1,1,1,1,1,1,
           0,1,1,1,1,2,2,2,
           0,1,2,2,2,2,2,2,
           0,1,2,3,3,3,3,3,
           0,1,2,3,3,3,3,4,
           0,1,2,3,4,4,4,4), nrow=8, byrow=TRUE)


  B = calculateLCS(A[1,],A[2,])
  expect_that( B, is_a("matrix"))
  expect_that( B, equals(result))
})

