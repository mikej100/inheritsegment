library(usethis)
library(tidyverse)
library(uuid)

#'
#' Make one crossover between chromesome a and b at location lock, List version
#'
#' chromosomes 1 and 2 are the parent chromosomes, chromosome 3 is the emerging
#' chilld chromosome. loc is the location 
#' @param ch_set list (
#'  a=input chromosome a,
#'  b=input chromosome b,
#'  c=result chromosome
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

#' Recombin a pair of chromosomes
#' 
#' Embedded in the function are parameters for the recombination:
#' exchange rate in often measured in centimorgans per megabyte of sequence,
#' here we use morgans per Mbase as the cx_rate, which is cM/100
#' 
#'
#' @param ch_pair chromosome a, chromosom b, chrom_number
#'
#' @return chromosome pair resulting from crossovers between chrom a and b.
#' @export
#'
#' @examples
recombine <- function (ch_pair, cx_rate = 0.01){  
  chrom_length <- last(ch_pair[[1]]$loc)
  generate_count <- 150
  
  
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
  if (sample( c(TRUE, FALSE),1)) {
    ch_0 <- list( a = ch_pair[[1]], b = ch_pair[[2]], c = make_ch( c(0), 1) ) 
  } else {
    ch_0 <- list( a = ch_pair[[2]], b = ch_pair[[1]], c = make_ch( c(0), 1) ) 
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
  chrom_a = tibble(
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
    chrom_a = mat_gamete_ch, 
    chrom_b = pat_gamete_ch,
    chrom_n = mat_chp$chrom_n
  )
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
# estimate average rate in Morgans (not cMorgans)
cross_rate_a <- 2 * genome_mbp / (genome_cm_female + genome_cm_male) /100
}

create_gen0_individual <- function (id, gender, population_id=1) {
  chr_number <- 3
   
  make_starting_ch_pair <- function(ch_num, id, ch_length = 130) {
    list(
      chrom_a =  make_ch(ch_length, id),
      chrom_b =  make_ch(ch_length, id),
      chrom_n = ch_num
      )
  }
  ch_set <- map(1:chr_number, \(ch_n) make_starting_ch_pair(ch_n, id))  
  
  individual <- list(
    id = id,
    gender = gender,
    pop_id = population_id,
    lineage <- id,
    ch_set = ch_set

  )
}

make_child <- function (mother, father) {
  
    make_child_ch_set <- function(mother, father) {
      map(seq_along(mother$ch_set),
                  ~ generate_zygote_chp(mother$ch_set[[.x]], father$ch_set[[.x]])
          )
    }
    child <- list (
      id = UUIDgenerate(),
      gender = sample(c("f", "m"), 1),
      pop_id = mother$pop_id,
      lineage = list(
        m_id = mother$id,
        p_id = father$id,
        m_lin =mother$lineage,
        p_lin = father$lineage
      ),
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
    ~ make_child(mothers[[1]], fathers[[1]])
  )               
}

make_generation <- function (generation) {
  sample_n <- length(generation) * 0.3
  mothers <- generation |> 
    keep(~ .x$gender =="f") |>
    sample(sample_n)
    
}
