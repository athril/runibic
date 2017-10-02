context("Calculating LCS")

test_that("Calculating Longest Common Subsequence (LCS): pairwiseLCS", {
  A = matrix(c(1,2,1,5,4,3,  2,1,3,2,1,4), nrow=2, byrow=TRUE)
  result=matrix(c(0,0,0,0,0,0,0,
                  0,0,1,1,1,1,1,
                  0,1,1,1,2,2,2,
                  0,1,2,2,2,3,3,
                  0,1,2,2,2,3,3,
                  0,1,2,2,2,3,4,
                  0,1,2,3,3,3,4), nrow=7, byrow=TRUE)


  B = pairwiseLCS(A[1,],A[2,])
  expect_that( B, is_a("matrix"))
  expect_that( B, equals(result))
})


test_that("Retrieving LCS: backtrackLCS", {
  A = matrix(c(1,2,1,5,4,3,  2,1,3,2,1,4), nrow=2, byrow=TRUE)
  result=c(1,2,1,4)
  C = backtrackLCS(A[1,],A[2,])
  expect_that(C, equals(result))
})
