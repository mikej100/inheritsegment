library(testthat)
library(purrr)

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("Make crossovers for given set of locations", {

 
  make_ch_0 <-function() {
    chrom_a <- make_ch( c(20, 40, 60, 80,100, 120, 130), 1)
    chrom_b <- make_ch( c(10, 30, 50, 70, 90, 110, 130), 1001)
    chrom_c <- make_ch( 0, 1)
    ch <- list(a=chrom_a, b=chrom_b, c=chrom_c)
  }

  # expected when xover at 65 and 105
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
test_that("Make recombination for given rate", {
  ch_length <- 130
  n <- 100
  cx_rate <- 0.01
  chrom_a = make_ch( ch_length, 1000)
  chrom_b = make_ch( ch_length, 2000)
  ch_pair <- list(
    chs = list(
      chrom_a,
      chrom_b),
    ch_n = 1
  )
  recombined <- map(  1:n, \(i) recombine(ch_pair, cx_rate = cx_rate ))
  cx_count <- map_int(recombined, \(r) nrow(r)) |>
    map_int( \(p) p-2)
  
  mean_cx_count <- mean(cx_count)
  poiss_exp <- cx_rate * ch_length
  sem <-  poiss_exp / n^0.5
  # Expect mean to be within four sigmas
  expect_equal(mean(cx_count), poiss_exp, tolerance = 4 * sem  )
  
})


test_that("Produces zygote chromosomes", {
  ch_length = 130
  m_chp <- list(
    chs = map(11:12, \(id) make_ch( ch_length, id )),
    ch_n = 1
  )
  p_chp <- list(
    chs = map(1:2, \(id) make_ch( ch_length, id )),
    ch_n = 1
  )
  xr_mat <- 0.08
  xr_pat <- 0.05
  
  n=3
  zygote_set <- map (1:n,
    \(i) generate_zygote_chp(m_chp, p_chp, xr_mat, xr_pat)
  )
  pat_ids <-zygote_set |> 
    map(\(z) z$chs[[2]]$id) |>
    reduce( \(cum, id) append(cum, id)) |>
    keep( \(id) !is.na(id))
  
  expect_equal(mean(pat_ids), 1.5, tolerance = 0.1 )  
  expect_equal(length(pat_ids), 130 * xr_pat * n, tolerance = 0.5)
})

zygote_set[[2]]$chs[[1]]$id
test_that("make population of generation 0 individuals", {
  females <- map(seq(1, length=12),
                 \(id) create_gen0_individual(id, gender = "f")
  )
  males <- map(seq(1001, length=10),
                 \(id) create_gen0_individual(id, gender = "m")
  )
  population <- c(females, males)
  
  expect_equal(population |>
                 keep( \(i) i$gender =="f") |>
                 length(),
               12
               )
  expect_equal(
    population |>
      keep (\(i) i$gender == "m") |>
      map_int(\(i) i$id) |>
      mean(),
    1005.5
  )
  mothers <<- sample(females, 5, replace=FALSE)
  fathers <<- sample(males, 5, replace=FALSE)
})

test_that("Produce child", {
  n_child <- 3
  child <- map(1, make_child(mothers[[1]], fathers[[1]]))
  
  children2 <- make_children(mothers[[1]], fathers[[1]], child_n = 3)
  children3 <- make_children(mothers[[1]], fathers[[1]])
  
  generation_1 <- map2(mothers, fathers,
                      \(m,f) make_children(m, f, child_n = 4)) |>
        reduce( ~ c(.x, .y))
  expect_equal(length(generation_1), 20)
  
  mothers2 <- generation_1 |>
    keep(~ .$gender == "f") |>
    sample(5, replace=FALSE)
  
  fathers2 <- generation_1 |>
    keep(~ .$gender == "m") |>
    sample(5, replace=FALSE)
  
  generation_2 <- map2(mothers2, fathers2,
                      \(m,f) make_children(m, f, child_n = 2)) |>
        reduce( ~ c(.x, .y))
  
  expect_equal(length(generation_2), 10)
  gen_2 <- as_tibble(transpose(generation_2))
})

test_that("Run generations",{
  pop_n <- 500
  gen_count <- 7
  growth_rate <- 1.1
  
  females <- map(seq(1, length=pop_n/2),
                 \(id) create_gen0_individual(id, gender = "f")
  )
  males <- map(seq(1001, length=pop_n/2),
                 \(id) create_gen0_individual(id, gender = "m")
  )
  gen_0 <- c(females, males)
  gens <- accumulate(1:gen_count,
                     \(gen, i) make_generation(gen, growth=growth_rate),
                     .init=gen_0 )  
  
  gen_size <- gens |>
    map_int(~length(.x))
  # Expect observed growth rate of about to match growth_rate parameter
  # Requires the generation size not to be randomised
  expect_equal(gen_size[3], 605)
  
  growth_rate_obs <-
    map2(gen_size, lag(gen_size),  ~ .x/.y ) |>
    tail(-2) |>
    unlist()
  expect_true (
    every( abs( growth_rate_obs - growth_rate ), ~ . <0.05)
  )
 
  gen_3 <- gens[[3]]
  gen_3_s <- sample (gens[[3]], 3)
  gen_3_ch_s <- get_chroms(gen_3_s )
  expect_equal(length(gen_3_ch_s), 6)
  
  # get list of all chromosomes in list of individuals
  get_chroms  <- function (gen) {
      chroms <- gen |>
      map( ~ .x$ch_set[[1]]$chs ) |>
      reduce( ~ c(.x, .y) )# |> 
  }
  
  # get all ids in all chromosom in list of chromosomes
  get_chroms_all_ids <- function (chroms) {
    chrom_ids <- chroms |>
    map(~ .x$id) |>
    reduce(~ c(.x, .y)) |>
    keep ( ~ !is.na(.x))
  }
  
  # get all chromosomes into one list, by generation
  gen_chroms <- gens |>
    map(~ get_chroms(.x))
  
  # get mean number of segments in a chromosome for each generation
  gens_id_n_mean <-  gen_chroms |>
    map_dbl( ~ .x|>
      map(~ .x$id) |>
      map_int( ~ length(.x) -1) |>
      mean()
    )
  gens_id_n_sd <-  gen_chroms |>
    map_dbl( ~ .x|>
      map(~ .x$id) |>
      map_int( ~ length(.x) -1) |>
      sd()
    )
  gens_id_n_growth   <-
      map2_dbl(gens_id_n_mean, lag(gens_id_n_mean), ~ .x/.y) 
  
  # get number of unique ids
  gens_unique_id_n <-  gen_chroms |>
    map( ~ .x|>
      map(~ .x$id)|>
      map( ~ unique(.x)) |>
      map_int( ~ length(.x) - 1 )
    )
  gens_unique_id_n_m <- map_dbl(gens_unique_id_n, ~ mean(.x)) 
  
  gens_unique_id_n_growth <-
    map2_dbl( gens_unique_id_n_m, lag(gens_unique_id_n_m), ~ .x/.y )
  get_ch_id_summary(gen_chroms[[1]])
  
  
  }
  
