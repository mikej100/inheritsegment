library(furrr)
library(tidyverse)
library(usethis)
library(uuid)

# Data structure =====
#, Data structures
#' in ABNF format
#' 
#'  generation = 1*individual
#'  individual = iid gender pop-id ch-set
#'  iid = INT | UUID ; individual identifier
#'  gender = "f" | "m"
#'  pop-id = INT  ; population id
#'  ch-set= 1*ch_pair ; chromosome set of the individual
#'  ch-pair = chs ch_n
#'  ch-n = INT   ; chromosome number
#'  chs = 2ch  ; chromosomes
#'  ch = loc id ; set of chromosome segments as lists of ids and locations
#'  ; length of id matches length of loc
#'  segment identifier, origin_id of the segment, NA marks end of chromosome
#'  id = 1*INT "NA"
#'  location on the chromosome of start of segment in Mbp last value end of chm
#'  loc = 0 1*REAL ; location on the chromosome of start of segment in Mbp
#'  
#'  
#'  lineage = 1*(iid-list)  ; first element is the first generation
#'  id-list = 1*(INT | UUID)  
#'  
#' Make one crossover between chromesome a and b at location lock. Dataframe ver
#'
#' Chromosome triplet
#'    cht$a is the chromsome being crossed into
#'    cht$b is the other chromosome, which is moved to position a before return
#'    cht$c is the new chromosome being built up
#' @param cht chromome triplet (a=parent 1, b=parent 1, c=child)
#'  @param location to crossover) 

#'
#' @return same strucure as param, with chrom 3 extended and chrom 1 and 2 
#' swapped 
#'
#' @examples
crossover <- function (cht, xloc) {
  # Logic of code
  #  loc
  #   all of c
  #   any of a > last of c && < xloc
  #   xloc
  #  
  #  id
  #   all of c less last c , which is NA
  #   a from last before c-last 
  #   to last. loc before xloc
  #   NA 
  
  # index of last location in $c
  i_c1 <- length(cht$c)    
  # index of last location in $a before end of $c
  i_a1 <- detect_index(cht$a$loc, ~ . > last(cht$c$loc))
  # index of last location in $a less than xloc
  i_a2 <- detect_index(cht$a$loc, ~ . >= xloc) -1 
  
  # all loc in $c
  n_loc1 <- cht$c$loc
  # any loc in a > last loc in c and < xloc
  n_loc2 <- head(cht$a$loc, i_a2  ) |>
    tail(i_a2 - i_a1 +1)
  
  # all ids from $c except terminating NA
  n_id1 <- head(cht$c$id, -1)
  # ids in $a from before end of c and up to xloc
  n_id2 <- head(cht$a$id, i_a2) |>
    tail(i_a2 - i_a1 +2  )
  
  new <- list(
    loc = c(n_loc1, n_loc2, xloc), 
    id =  c(n_id1, n_id2, NA)
  )
  
 ch_return <- list (
    a = cht$b,
    b = cht$a,
    c = new
  )
}
#' Recombine a pair of chromosomes
#' 
#' Embedded in the function are parameters for the recombination:
#' exchange rate in often measured in centimorgans per megabyte of sequence,
#' here we use morgans per Mbase as the cx_rate, which is cM/100
#' 
#'
#' @param ch_pair list(chrom a, chromo b), chromosome number
#'
#' @return chromosome pair resulting from crossovers between chrom a and b.
#' @export
#'
#' @examples
recombine <- function (ch_pair, cx_rate = 0.01){  
  chrom_length <- last(ch_pair$chs[[1]]$loc)
  generate_count <- 50
  
  
  locs_1 <- rexp(generate_count, cx_rate) |>
    accumulate( \(len, loc) len + loc ) 
  
  # Will the generated locations be long enough?
  # If this stop is thrown, then increase the generate_count.
  # map_dbl( c(0.5, 0.1, .1e-3, 1e-12, 1e-15),\(p) stats::qgamma( p, shape = 150, rate = .01 ))
  min_length <- stats::qgamma( 1e-15, shape = generate_count, rate = cx_rate)
  if (min_length < chrom_length) {
    stop("Warning: generating crossover location list may not be long enough")
  }
  if (last(locs_1) <= chrom_length)  {
    stop("Error generating crossover location list, too short")
  }
  
  locs <- head_while(locs_1, \(len) len <= chrom_length) |>
    append(chrom_length)
  
  # randomise choice of which chromosome used as "a"
  if (sample( c(TRUE, FALSE), 1)) {
    ch_0 <- list( a = ch_pair$chs[[1]],
                  b = ch_pair$chs[[2]], 
                  c = make_ch( c(0), 1) ) 
  } else {
    ch_0 <- list( a = ch_pair$chs[[2]],
                  b = ch_pair$chs[[1]],
                  c = make_ch( c(0), 1) ) 
  }
   
  reduce( locs, 
                        \(ch, xl) crossover (ch, xl),
                        .init=ch_0
                        )$c
}
    

#' Set up a chromosome with given set of segments with sequential ids
#' 
#' Can be used to set up multiple segments or whole chromosome with single id.
#' Special case to create an empty chromosome locs <- 0
#' @param locs List of segment end positions, including end of chromosome
#' @param start_id Start number for the segment ids
#'
#' @return Tibble specifying the segments
#'
#' @examples
make_ch <- function (locs, start_id) {
  if (locs[1] == 0){
    locsm <- 0
  } else {
    locsm <- c(0, locs)
  }
  ch = list(
    loc = locsm,
    id = c(seq(start_id, length= length(locsm) - 1), NA)
  )
}

#' Generate zygote chromosome pair from maternal and paternal chromosome pairs
#' i.e. meiosis with crossover, selection and fusion.
#'
#' @param mat_ch maternal chromosome pair
#' @param pat_ch paternal chromosome pair
#' @param mat_xr maternal crossover rate
#' @param pat_xr paternal crossover rate
#'
#' @returno 
#' @export zygote chromosome pair
#'
#' @examples
generate_zygote_chp <- function (mat_chp, pat_chp, mat_xr=0.01, pat_xr=0.01 ){
  mat_gamete_ch <- recombine(mat_chp, mat_xr)
  pat_gamete_ch<- recombine(pat_chp, pat_xr)
  
  zygote_chp <- list(
    chs = list(
      mat_gamete_ch, 
      pat_gamete_ch
      ),
    ch_n = mat_chp$ch_n
  )
}   

#   make_chrome <- function (type) {
#   chrom_a |>
#     head_while(\(p) p$loc <= c_loc)
#   # chromosome is made up of segments from ancestors
#   # has fixed length
#   # segment has start and end location, and id of root ancestor
#   
#   
#   # crossover rate
#   # female genome is 4460 cM, and male is 2590 cm (wikipedia)
#   # Take human genome 3.1 e9 long, in 24 chromosomes
#   genome_cm_female <- 4460
#   genome_cm_male <- 2590
#   genome_mbp <- 3.1e3
#   chrom_count <- 24
#   
#   chrom_size <- genome_mbp / 24
#   # estimate average rate in Morgans (not cMorgans)
#   cross_rate_a <- 2 * genome_mbp / (genome_cm_female + genome_cm_male) /100
#   }

create_gen0_individual <- function (id, gender, ch_count=1, population_id=1) {
   
  make_starting_ch_pair <- function(ch_num, id, ch_length = 130) {
    list(
      chs = list (
        make_ch(ch_length, id),
        make_ch(ch_length, id)
      ),
      ch_n = ch_num
      )
  }
  ch_set <- map(1:ch_count, \(ch_n) make_starting_ch_pair(ch_n, id))  
  
  individual <- list(
    iid = id,
    gender = gender,
    pop_id = population_id,
    lineage = list(id),
    ch_set = ch_set

  )
}

make_child <- function (mother, father) {
  
    make_child_ch_set <- function(mother, father) {
      map(seq_along(mother$ch_set),
                  ~ generate_zygote_chp(mother$ch_set[[.x]], father$ch_set[[.x]])
          )
    }
    iid <- UUIDgenerate()
    
    lineage <- map2( mother$lineage, father$lineage, ~  c(.x, .y) ) |>
      c(iid)

    child <- list (
      gender = sample(c("f", "m"), 1),
      iid = iid,
      pop_id = mother$pop_id,
      lineage = lineage,
      ch_set = make_child_ch_set(mother, father)
    )
}

make_children <- function(mother, father, child_n = NA, child_mean=2) {
  # Randomise the number of children, if not specified.
  if (is.na(child_n)){
    child_n <- rpois(1, child_mean)
  }
  
  children <- map(
    seq(1, length=child_n),
    ~ make_child(mother, father)
  )               
}
furrr_generated <- function(mothers, fathers, child_m) {
  # Split generation for furrr
  future::plan(multisession, workers=6)
  generated_a <- future_map2(mothers, fathers, 
                      \(m,f) make_children(m, f, child_n = child_m),
                      .options = furrr_options(seed=TRUE)
                      )
  generated <- generated_a |>
    reduce( ~ c(.x, .y))
}

make_generation <- function (parent_gen, growth=1.1, limit_m=NA) {
  child_m <- 6  #mean number of children in family
  
  # estimate quantity breeding pairs needed based on target growth, with margin
  par_n <- length(parent_gen)
  if (is.na(limit_m)) {
    sample_n_target <- par_n * growth / child_m * 1.1
  } else {
    # TODO calculate samples size when limit is defined
    sample_n_target <- par_n
  }
  # Use code like this if want to randomise breeding pairs
  # sample_n <- rpois(1, length(parent_gen)/2*breed_rate)
  
  females <- parent_gen |> 
    keep(~ .x$gender =="f")
  males <- parent_gen |>
    keep(~ .x$gender =="m")
  
  sample_n_available <- min(length(females), length(males))
  sample_n <- min(sample_n_available, sample_n_target)
  
  mothers <- sample(females, sample_n , replace = TRUE )
  fathers <- sample(males, sample_n , replace = TRUE )   
  
  furrr <- T
  if (furrr) {
    generated <- furrr_generated(mothers, fathers, child_m)
  } else {
  generated <- map2(mothers, fathers,
                      \(m,f) make_children(m, f, child_n = child_m)) |>
    reduce( ~ c(.x, .y))
  }
  if (is.na(limit_m)) {
    limit <- length(parent_gen) * growth 
  }
  
  gen_next <- sample(generated, limit) 
  # use this to randomise generation size
  # gen_next <- sample(generated, rpois(1, limit)) 
}
