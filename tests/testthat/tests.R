context(
  "Compare results with  Marcotte, Denis. 1996. 'Fast Variogram Computation with FFT.' Computers & Geosciences 22 (10): 1175–86. doi:10.1016/S0098-3004(96)00026-X."
)

## set up input matrices from paper
m1 = raster(matrix(c(3, 6, 5, 7, 2, 2, 4, NA, 0), ncol = 3, byrow = T))
m2 = raster(matrix(c(10, NA, 5, NA, 8, 7, 5, 9, 11), ncol = 3, byrow = T))

ac = acorr(m1, padlongitude = T, verbose = T)

test_that("Confirm calculation of number of observations==nh11 on top of page 1179",
          {
            nh11=matrix(
              c(1,1,2,1,1,
                2,3,5,3,2,
                3,4,8,4,3,
                2,3,5,3,2,
                1,1,2,1,1),
              nrow=5,byrow=T)

          nh11_test=10^as.matrix(ac[["nobs"]])
          expect_equal(nh11,nh11_test)
          })

test_that("Confirm results from gh11 on top of page 1179", 
          {
            gh11=matrix(c(
               0.000,  0.000, -2.000,  0.000,  0.000,
              -2.000,  1.111,  0.400,  1.222,  2.250,
              -1.222, -0.375,  4.734, -0.375, -1.222,
               2.250,  1.222,  0.400,  1.111, -2.000,
               0.000,  0.000, -2.000,  0.000,  0.000),
              nrow=5,byrow=T)
            gh11_test = round(as.matrix(ac[["acor"]]) / 10, 3)
            expect_equal(gh11, gh11_test)
})
