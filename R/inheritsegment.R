library(furrr)
library(logger)
library(lubridate)
library(tictoc)
library(tidyverse)
library(usethis)
library(uuid)

log_info("Logging inheritsegment")
# Data structure =====
#, Data structures
#' in ABNF format
#' 
#'  generation = 1*individual ; same as a population
#'  individual = iid gender pop-id ch-set
#'  iid = INT | UUID ; individual identifier
#'  gender = "f" | "m"
#'  subpop-id = INT  ; sub-population id
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
#'  lineage = 1*(iid-list)  ; list for each generation, starting with oldest
#'  iid-list = 1*(INT | UUID)   ; is of individual id from that generation

#'  
#' Make one crossover between chromesome a and b at location loc
#'
#' @param cht chromosome triplet 
#'    cht$a is the chromsome being crossed into
#'    cht$b is the other chromosome, which is moved to position a before return
#'    cht$c is the new chromosome being built up
#'  @param xloc location to crossover at
#'  
#' Logic of code to build up chromosome c
#'  loc:
#'   all of c
#'   any of a > last-of-c && < xloc
#'   xloc
#'  
#'  id:
#'   all of c less last c , which is NA
#'   a from last before last-of-c to last loc before xloc
#'   NA 
#'
#' @return same cht structure as param, with chrom c extended and chrom a and b
#' swapped
#'
#' @examples
crossover <- function (cht, xloc) {
  
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
  generate_count <- 5000 * cx_rate
  
  # generate random locations for crossover
  locs_1 <- rexp(generate_count, cx_rate) |>
    accumulate( \(len, loc) len + loc ) 
  
  # Will the generated locations be long enough?
  # TODO replace all this checking with a loop to generate more if needed.
  # If this stop is thrown, then increase the generate_count.
  # Use gamma function to assess whether probability of failure is very low
  # map_dbl( c(0.5, 0.1, .1e-3, 1e-12, 1e-15),\(p) stats::qgamma( p, shape = 150, rate = .01 ))
  min_length <- stats::qgamma( 1e-15, shape = generate_count, rate = cx_rate)
  if (min_length < chrom_length) {
    # warn if there is perceptible probability of failure
    stop("Warning: generating crossover location list may not be long enough")
  }
  if (last(locs_1) <= chrom_length)  {
    # It happened, not enough locs, need to increase
    stop("Error generating crossover location list, generate_count too small")
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
#' @return list specifying the segments
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
#' @param xr crossover rates vector of maternal and paternal rates. Crossover 
#' per Mbp default 0.0144 and 0.083 (human average)
#'
#' @return zygote chromosome pair
#'
#' @examples
#' 
generate_zygote_chp <- function (mat_chp, pat_chp, xr=c(0.0144, 0.0083)){
  mat_gamete_ch <- recombine(mat_chp, xr[1])
  pat_gamete_ch <- recombine(pat_chp, xr[2])
  
  zygote_chp <- list(
    chs = list(
      mat_gamete_ch, 
      pat_gamete_ch
      ),
    ch_n = mat_chp$ch_n
  )
}   

#   
#   # crossover rate
#   # female genome is 4460 cM, and male is 2590 cm (wikipedia)
#   # Take human genome 3.1 e9 long, in 23 chromosomes
#   genome_cm_female <- 4460
#   genome_cm_male <- 2590
#   genome_mbp <- 3.1e3
#   chrom_count <- 23
#
#      
#   chrom_size <- genome_mbp / chrom_count = 135
#   maternal and paternal crossover rate average 
#   4460/3100 = 0.0144 paternal 2590/3100 = 0.00835
#   # estimate average rate in Morgans (not cMorgans)
#   cross_rate_m <- genome_mbp / (genome_cm_female /100
#   }

#' Create first generation individual, ancestors of following generations
#'
#' @param id integer, id for segments of all chromosomes, also is individual id
#' @param gender "f" or "m" 
#' @param ch_count integer, haploid count of chromosomes individual to have. 
#' default=2 
#' @param subpop_id integer, id to specify a sub-population in the generation. 
#' default=1k
#'
#' @return list with structure of an individual
#'
#' @examples
create_gen1_individual <- function (id, gender, ch_count=2, subpop_id=1) {
   
  make_starting_ch_pair <- function(ch_num, id, ch_length = 135) {
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
    subpop_id = subpop_id,
    lineage = list(id),
    ch_set = ch_set

  )
}

#' Generate child from two parents
#'
#' Create new individual with chromosome set derived from mother and father
#' meiosis, chromosome selection to haploid, and fusion to diploid. Randomise
#' gender and assign UUID to individual id.
#'
#' @param mother individual which is the mother 
#' @param father individual which is the father 
#'
#' @return new individual with set of chromosomes derived from parents
#'
#' @examples
make_child <- function (mother, father) {
  
    make_child_ch_set <- function(mother, father) {
      map(seq_along(mother$ch_set),
                  ~ generate_zygote_chp(mother$ch_set[[.x]], father$ch_set[[.x]])
          )
    }
    iid <- UUIDgenerate()
    
     
    lineage <- c(mother$iid, father$iid)
#    lineage <- map2( mother$lineage, father$lineage, ~  c(.x, .y) ) |>
#      c(iid)

    child <- list (
      gender = sample(c("f", "m"), 1),
      iid = iid,
      subpop_id = mother$subpop_id,     # Simply assign mother's subpopulation id
      lineage = lineage,
      ch_set = make_child_ch_set(mother, father)
    )
}

#' Generate a number of children from pair of parents
#'
#' @param mother individual which is the mother 
#' @param father individual which is the father 
#' @param child_n (optional) number of children for the couple
#' @param child_mean mean number of children for the couple, used to randomised
#' the number of children. 
#'
#' @return list of new individuals with set of chromosomes derived from parents
#'
#' @examples
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

# furrr code isolated in its own function to minimise objects in environment 
furrr_generated <- function(mothers, fathers, child_m) {
  # Split generation for furrr
  future_workers <- 6 
  log_info("Parallel processing future_workers: {future_workers}")
  future::plan(multisession, workers=future_workers)
  generated_a <- future_map2(mothers, fathers, 
                      \(m,f) make_children(m, f, child_n = child_m),
                      .options = furrr_options(seed=TRUE)
                      )
  generated <- generated_a |>
    reduce( ~ c(.x, .y))
}

#' Create a new generation of individuals from a given generation.
#'
#' @param parent_gen parent generation
#' @param growth factor by which population size will increase each generation
#' @param limit_m limit of size of generation (not implemented) 
#' @param parallel set TRUE to run in parallel (for large populations), default
#' is FALSE (faster for small tests)
#'
#' @return new generation
#'
#' @examples
make_generation <- function (parent_gen, growth=1.1, limit_m=NA,
                             parallel=FALSE) {
  child_m <- 3  #mean number of children in family
  log_info("Making generation") 
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
  
  if (parallel) {
    generated <- furrr_generated(mothers, fathers, child_m)
  } else {
  generated <- map2(mothers, fathers,
                      \(m,f) make_children(m, f, child_n = child_m)) |>
    reduce( ~ c(.x, .y))
  }
  if (is.na(limit_m)) {
    limit <- length(parent_gen) * growth 
  }
  
  # randomise order and selection for next generation
  gen_next <- sample(generated, limit) 
  # use this to randomise generation size
  # gen_next <- sample(generated, rpois(1, limit)) 
}

  generate <- function(  pop_n,  gen_count, ch_count = 23, growth_rate = 1.0) {
 
  gen_1 <- c(
        map(seq(1, length=pop_n/2),
                 \(id) create_gen1_individual(
                   id, gender = "f", ch_count = ch_count
                 )
        ),
        map(seq(500001, length=pop_n/2),
                 \(id) create_gen1_individual(
                   id, gender = "m", ch_count = ch_count)
            )
  )
  parallel_switch <- pop_n * ch_count >= 3000
  tic()
  log_info("Generation parameters: 
            starting population: {pop_n}
            generation count: {gen_count}
            growth rate: {growth_rate}
            chromsome count: {ch_count}
            parallel switch: {parallel_switch}
           " )
  gens <- accumulate(1:gen_count, .init=gen_1, 
                     \(gen, i) make_generation(
                                    gen,
                                    growth=growth_rate,
                                    parallel = parallel_switch
                                    )
                     )
  timer <- toc()
  log_info("Generation timer: {timer$callback_msg}") 
  
  return(gens)
  }
  
  write_data <- function(data, filename) {
    full_filename <- str_c(filename,".RData")
    filepath <- file.path("..","not_in_repo","data","generations",full_filename)
    save(data, file=filepath)
  }
  load_data <- function(filename) {
    full_filename <- str_c(filename,".RData")
    filepath <- file.path("..","not_in_repo","data","generations",full_filename)
    result <- load(filepath)
    
  }
# Extract timing intervals from logfile  
log_timing_analysis <- function(text){
  
  datetimes <- function (text) {
    str_extract_all(text,"\\d{4}-\\d{2}-\\d{2}.{9}") |>
      unlist()
  }   
  dt <- datetimes(text)
}
  