library(testthat)

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

make_ch <- function (locs, start_id) {
  chrom_a = tibble(
    loc = locs,
    id = c(seq(start_id, length= length(locs)-1), NA)
  )
}

test_that("Recombine along a chromosome", {

 
  make_ch_0 <-function() {
    chrom_a <- make_ch( c(0, 20, 40, 60, 80,100, 120, 130), 1)
    chrom_b <- make_ch( c(0, 10, 30, 50, 70, 90, 110, 130), 1001)
    chrom_c <- make_ch( c(0), 1)
    ch <- list(a=chrom_a, b=chrom_b, c=chrom_c)
  }

  # expected if xover at 65 and 105
  expected_1 <- tibble(
    loc = c(0, 20, 40, 60,   65,   70,   90, 105, 120, 130),
    id =  c(1,  2,  3,  4, 1004, 1005, 1006,   6,   7, NA)
  )
  
  ch_0 <- make_ch_0()
  recombined <- reduce( c(65, 105, 130), 
                        \(ch, xl) crossover (ch, xl),
                        .init=ch_0
                        )$c
  expect_equal(recombined, expected_1)
  # res <- accumulate(c(65, 95, 130), \(ch, xl) crossover(ch, xl), .init=ch)
})



# For testing list versions of recombination, not used ----
# make_ch_l <-function() {
#   chrom_a = list(
#     loc = c(0, 20, 40, 60, 80, 130),
#     id = c(seq(1, length=5), NA)
#   ) |>
#     transpose()
#   chrom_b = list(
#     loc = c(0, 10, 30, 50, 70, 90, 130),
#     id = c(seq(6, length=6), NA)
#   ) |>
#     transpose()
#   chrom_c <- list(
#     loc = 0.0,
#     id = NA
#   ) |>
#     transpose()
#   ch <- list(a=chrom_a, b=chrom_b, c=chrom_c)
# })