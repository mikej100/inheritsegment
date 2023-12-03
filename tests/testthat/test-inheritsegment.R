library(testthat)
library(dplyr)
library(purrr)

test_that("Make crossovers for given set of locations", {
 
  make_cht_1 <- function() {
    chrom_a <- make_ch( c(20, 40, 60, 80,100, 120, 130), 1)
    chrom_b <- make_ch( c(10, 30, 50, 70, 90, 110, 130), 1001)
    chrom_c <- make_ch( 0, 1)
    cht <- list(a=chrom_a, b=chrom_b, c=chrom_c)
  }

  # expected when xover at 65 and 105
  expected_1 <- list(
    loc = c(0, 20, 40, 60,   65,   70,   90, 105, 120, 130),
    id =  c(1,  2,  3,  4, 1004, 1005, 1006,   6,   7, NA)
  )
  
  cht_1 <- make_cht_1()
  recombined <- reduce( c(65, 105, 130), 
                        \(cht, xl) crossover (cht, xl),
                        .init=cht_1
                        )$c
  expect_equal(recombined, expected_1)
})

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
  cx_count <- map_int(recombined, ~ length(.$loc)) |>
    map_int( \(p) p-2)
  
  mean_cx_count <- mean(cx_count)
  poiss_exp <- cx_rate * ch_length
  sem <-  poiss_exp / n^0.5
  # Expect mean to be within four sigmas
  expect_true( abs(mean(cx_count) - poiss_exp) < 4 * sem  )
  
})

test_that("Produce zygote chromosome pairs", {
  mother <- map(1:2, ~ create_gen1_individual(.x, "f") )
  father<- map(1:2, ~ create_gen1_individual(1000 * .x + 1 , "m") )
  
  # create first generation
  m1_chp <- mother[[1]]$ch_set[[1]]
  p1_chp <- father[[1]]$ch_set[[1]]
  xr <- c(0.1, 0.05)
  zygote1_chp <-  generate_zygote_chp(m1_chp, p1_chp, xr)
  
  # check first generation 
  pat1_id <- zygote1_chp$chs[[2]]$id |> head(-1)
  expect_true( every(pat1_id, ~ .x == 1001))
  
  mat1_id <- zygote1_chp$chs[[1]]$id |> head(-1)
  expect_gte (length(mat1_id), length(pat1_id) )
  
  # create second generation
  m2_chp <- mother[[2]]$ch_set[[1]]
  p2_chp <- father[[2]]$ch_set[[1]]
  zygote2_chp <-  generate_zygote_chp(m2_chp, p2_chp, xr)
  zygote3_chp <-  generate_zygote_chp(zygote1_chp, zygote2_chp, xr)
  
  # check second generation
  pat3_id <- zygote3_chp$chs[[2]]$id |> head(-1)
  expect_true( every(pat3_id, ~ .x %in% c(2, 2001)))
}) 

test_that("create individual of generation 1",{
  ind1 <- create_gen1_individual(2, gender = "f")
  expect_equal( length(ind1$ch_set), 2)
  
  ind2 <- create_gen1_individual(2, gender = "m", ch_count = 5)
  expect_equal( length(ind2$ch_set), 5)
})

test_that("make population of generation 1 individuals", {
  females <- map(seq(1, length=12),
                 \(id) create_gen1_individual(id, gender = "f")
  )
  males <- map(seq(1001, length=10),
                 \(id) create_gen1_individual(id, gender = "m")
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
      map_int(\(i) i$iid) |>
      mean(),
    1005.5
  )
})

test_that("Produce children and generations", {
  mothers <- map(seq(1, length=5),
                 \(id) create_gen1_individual(id, gender = "f")
  )
  fathers <- map(seq(1001, length=5),
                 \(id) create_gen1_individual(id, gender = "m")
  )
  
  # check child has expected ids
  child <- make_child(mothers[[2]], fathers[[3]])
  child_mat_id <- child$ch_set[[2]]$chs[[1]]$id |> head(-1)
  child_pat_id <- child$ch_set[[2]]$chs[[2]]$id |> head(-1)
  expect_true( every(child_mat_id, ~ . == 2))
  expect_true( every(child_pat_id, ~ . == 1003))
  
  # check lineage as expected
  expect_equal(child$lineage[[1]], c(2,1003))
  # expect uuid
  expect_equal(str_length(child$lineage[[2]]) , 36 )
  
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
  
  expect_equal(length(generation_1[[3]]$lineage), 2)
  expect_equal(length(generation_2[[3]]$lineage), 3)
})

test_that("Run generations",{
  # set to 100 for these tests
  pop_n <- 100
  # set gen_count to at least 3 for the tests
  gen_count <- 3
  growth_rate <- 1.1
  
  # females <- map(seq(1, length=pop_n/2),
  #                \(id) create_gen0_individual(id, gender = "f")
  # )
  # males <- map(seq(1001, length=pop_n/2),
  #                \(id) create_gen0_individual(id, gender = "m")
  # )
  
  #create males and females in one step to limit object size in environment for 
  # furrr
  gen_1 <- c(
        map(seq(1, length=pop_n/2),
                 \(id) create_gen1_individual(id, gender = "f")
        ),
        map(seq(1001, length=pop_n/2),
                 \(id) create_gen1_individual(id, gender = "m")
            )
  )
  
  # Generate multiple generations ====================
  gens <- accumulate(1:gen_count,
                     \(gen, i) make_generation(gen, growth=growth_rate),
                     .init=gen_1 )  
  
  # Expect observed growth rate of about to match growth_rate parameter
  # Requires the generation size not to be randomised
  gen_size <- gens |>
    map_int(~length(.x))
  expect_equal(gen_size[3], 121)
  
  growth_rate_obs <-
    map2(gen_size, lag(gen_size),  ~ .x/.y ) |>
    tail(-2) |>
    unlist()
  expect_true (
    every( abs( growth_rate_obs - growth_rate ), ~ . <0.05)
  )
  
  # Look at non-unique iids in lineage, result of related ancestry)
  sample_n <- 50
  non_unique_ancestor_rate <- sample(last(gens), sample_n )  |> 
    map(~ .$lineage[[1]])  |>
    map(~ unique(.)) |>
    map_int(~ length(.))|>
    keep( ~ . != 2^gen_count) |>
    length()/sample_n
  expect_lt( abs(non_unique_ancestor_rate - 0.28), .3) 
  
  
  # get list of all chromosomes 1 in list of individuals
  get_chroms  <- function (gen) {
      chroms <- gen |>
      map( ~ .x$ch_set[[1]]$chs ) |>
      reduce( ~ c(.x, .y) )# |> 
  }
 
  # test can extract expected number of chromosomes 
  gen_3_s <- sample (gens[[3]], 3)
  gen_3_ch_s <- get_chroms(gen_3_s )
  expect_equal(length(gen_3_ch_s), 6)
  
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
  # expect mean number of chromosome segments to grow continuously through gens
  # list of increments in ids
  id_inc   <-
      map2_dbl(gens_id_n_mean, lag(gens_id_n_mean), ~ .x-.y) |>
    tail( -1)
  
  id_inc_expected <- 1.2
  id_inc_deviation <- map_dbl(id_inc, ~ abs(.x - id_inc_expected))
  })

test_that("Model the chances of sharing any DNA with someone n-generations back", {
  pop_n <- 1e2
  gen_count <- 4
  growth_rate <- 1.
  
  gen_1 <- c(
        map(seq(1, length=pop_n/2),
                 \(id) create_gen0_individual(id, gender = "f")
        ),
        map(seq(2001, length=pop_n/2),
                 \(id) create_gen0_individual(id, gender = "m")
            )
  )
  # Generate multiple generations ====================
  parallel_switch <- pop_n >= 1000
  gens <- accumulate(1:gen_count, .init=gen_1, 
                     \(gen, i) make_generation(
                                    gen,
                                    growth=growth_rate,
                                    parallel = parallel_switch
                                    )
                     )
  # get all chromosomes into one list, by generation
  gen_chroms <- gens |>
    map(~ get_chroms(.x))
  
  
  
  # get number of unique ids per chromosome set for each individual, by generation
  gen_unique_id_n <-  gen_chroms |>
    map( ~ .x|>
      map(~ head(.x$id, -1)) |>
      map( ~ unique(.x)) |>
      map_int( ~ length(.x)  )
      )
  gen_chroms[[2]] |> map(~ head(.x$id, -1)) |> sample(10) 
  sample(gen_unique_id_n[[3]],20)
  gens_unique_id_n_m <- map_dbl(gens_unique_id_n, ~ mean(.x)) 
  
  gens_unique_id_n_growth <-
    map2_dbl( gens_unique_id_n_m, lag(gens_unique_id_n_m), ~ .x/.y )
  
})
