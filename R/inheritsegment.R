library(usethis)
library(tidyverse)

#'
#' Make one crossover between chromesome a and b at location lock, List version
#'
#' chromosomes 1 and 2 are the parent chromosomes, chromosome 3 is the emerging
#' chilld chromosome. loc is the location 
#' @param ch_set list (
#'  input chromosome 1,
#'  input chromosome 2,
#'  result chromosome
#'  @param location to crossover) 

#'
#' @return same strucure as param, with chrom 3 extended and chrom 1 and 2 
#' swapped 
#'
#' @examples
crossover_l <- function (ch, xloc) {

 segment_head <- ch$a |>
    tail_while( \(cha) cha$loc >= last(ch$c)$loc) 
  segment <- segment_head |>
    head_while( \(sh) sh$loc <= xloc )
  
  last <- list(list(loc=xloc, id = NA))
  c_return1 <- append(ch$c, segment)
  c_return2 <- append(c_return1, last)
  
  c_return3 <- append(ch$c, c(segment, last))
  #last_id <- detect(ch$a, \(cha) cha$loc < xloc, .dir="backward")
  ch_return <- list (
    a = ch$b,
    b = ch$a,
    c = append(ch$c, c(segment, last))
  )
  return (ch_return) 
}
#'
#' Make one crossover between chromesome a and b at location lock. Dataframe ver
#'
#' chromosomes 1 and 2 are the parent chromosomes, chromosome 3 is the emerging
#' chilld chromosome. loc is the location 
#' @param ch_set list (
#'  input chromosome 1,
#'  input chromosome 2,
#'  result chromosome
#'  @param location to crossover) 

#'
#' @return same strucure as param, with chrom 3 extended and chrom 1 and 2 
#' swapped 
#'
#' @examples
crossover <- function (ch, xloc) {
  start_loc <- last(ch$c$loc)  
  leader <- head(ch$c, -1)
  
  start_ind <- detect_index(ch$a$loc, \(l) l > start_loc)  -1 
  start_id <- ch$a$id[[start_ind]]
  first <- tibble(loc=start_loc, id=start_id)
  
  segment <- ch$a |>
    filter(loc > start_loc & loc < xloc) 
  
  last <- tibble(loc=xloc, id=NA)
  
  ch_return <- list (
    a = ch$b,
    b = ch$a,
    c = rbind(leader, first, segment, last)
  )
}

recombine <- function (chrom_a, chrom_b){  
  chrom_length <- 10
  cx_rate <- 0.8
  generate_length <- 300
  
  locs_1 <- rexp(generate_length, cx_rate) |>
    accumulate( \(len, loc) len + loc ) 
  # TODO use the gamma approximation to sum of exponentials to specify
  # a safe number of locations
  if (last(locs_1) <= chrom_length)  {
    stop("Error generating crossover location list")
    # use stats::qgamma( 1 - 1e-15, shape = gen_len, rate = cx_rate)
    # to check shortest length for a trillion tries
  }
}
  
  

make_chrome <- function (type) {
chrom_a |>
  head_while(\(p) p$loc <= c_loc)
# chromosome is made up of segments from ancestors
# has fixed length
# segment has start and end location, and id of root ancestor


# crossover rate
# female genome is 4460 cM, and male is 2590 cm (wikipedia)
# Take human genome 3.1 e9 long, in 24 chromosomes
genome_cm_female <- 4460
genome_cm_male <- 2590
genome_mbp <- 3.1e3
chrom_count <- 24

chrom_size <- genome_mbp / 24
cross_rate_a <- 2 * genome_mbp / (genome_cm_female + genome_cm_male)
}