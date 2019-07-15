#!/usr/bin/env R
# TCRMaxSAT
# Copyright (C) 2019    Rahul Dhodapkar

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

library(hashmap)
library(Matrix)
library(dplyr)

strreverse <- function(x){ strsplit(x, NULL) %>% lapply(rev) %>% sapply(paste, collapse="") }

residue_residue_rt_energies <- 
  as.matrix(read.csv('./config/residue_residue_rt_energies.csv', row.names = 1))
class(residue_residue_rt_energies) <- "numeric"
residue_residue_rt_energies <- forceSymmetric(residue_residue_rt_energies)

respairs <- c()
energies <- c()
for(r1 in rownames(residue_residue_rt_energies)) {
  for(r2 in colnames(residue_residue_rt_energies)) {
    respairs <- c(respairs, paste0(r1, r2))
    energies <- c(energies, residue_residue_rt_energies[r1,r2])
  }
}
respair2energy <- hashmap(respairs, energies)
rm(respairs, energies)

aa_info <- 
  read.csv('./config/aa_info.csv', header=T)

olc2tlc <- hashmap(as.character(aa_info$One.Letter.Code), as.character(aa_info$Three.Letter.Code))
ucname2tlc <- hashmap(toupper(as.character(aa_info$Name)), as.character(aa_info$Three.Letter.Code))

#' Affinity score between amino acid sequences.
#'
#' \code{affinity} returns an affinity score between affinity acid sequences
#'
#' @param seq.1 character vector; first amino acid sequence
#' @param seq.2 character vector; second amino acid sequence
#' @param code.type string, one of 'olc', 'tlc', or 'name', defining
#'   in what form the amino acid codes will be in \code{sequence}
#' @param binding.width Integer. the number of amino acids in the modeled CDR3-epitope contact
#' @return an affinity score that is determined by the residue-residue
#'   energies of amino acids in the sequence
#' @examples
#' affinity(c('I', 'L', 'G', 'G', 'G', 'G', 'G'), 
#'          c('G', 'G', 'G', 'G', 'I', 'L', 'C', 'G'))
#' affinity('ILGGGCCCQL', 'GGGGTYFPGLI')
#' \dontrun{
#' affinity('ILGG', 'GGGG')
#' }
affinity <- function(seq.1, seq.2, code.type='olc', binding.width = 3) {
  # check input types
  if (is.character(seq.1) && length(seq.1) == 1 && code.type =='olc') {
    seq.1 <- unlist(strsplit(seq.1, ''));
    seq.2 <- unlist(strsplit(seq.2, ''));
  }
  
  # convert to tlc sequence
  if (code.type == 'olc') {
    seq.1 <- olc2tlc[[seq.1]];
    seq.2 <- olc2tlc[[seq.2]];
  } else if (code.type == 'tlc') {
    # Do nothing
  } else if (code.type == 'name') {
    seq.1 <- ucname2tlc[[toupper(seq.1)]];
    seq.2 <- ucname2tlc[[toupper(seq.2)]];
  }
  
  min.offset.1 <- 1;
  min.offset.2 <- 1;
  min.rt.energy <- 0;
  
  for (i.1 in seq(1, length(seq.1) - binding.width)) {
    for (i.2 in seq(1, length(seq.2) - binding.width)) {
      nrg <- sum(respair2energy[[ 
          paste0(seq.1[i.1:(i.1+binding.width)],
                 seq.2[i.2:(i.2+binding.width)])
        ]])

      if (nrg < min.rt.energy) {
        min.offset.1 <- i.1;
        min.offset.2 <- i.2;
        min.rt.energy <- nrg;
      }
    }
  }
  
  return(list(
    min.offset.1 = min.offset.1,
    min.offset.2 = min.offset.2,
    min.rt.energy = min.rt.energy
  ));
}


#' Generate random amino acid sequence of defined length
#'
#' \code{random_aa} returns a random amino acid sequence
#'
#' @param n Integer length of the amino acid sequence desired
#' @param code.type string, one of 'olc', 'tlc', or 'name', defining
#'   in what form the amino acid codes will be
#' @param binding.width Integer. the number of amino acids in the modeled CDR3-epitope contact
#' @return a character vector of an amino acid sequence
#' @examples
#' random_aa(5)
random_aa <- function(n = 7, code.type = 'olc') {
  if (code.type == 'olc') {
    char_bank <- as.character(aa_info$One.Letter.Code)
  } else if (code.type == 'tlc') {
    char_bank <- as.character(aa_info$Three.Letter.Code)
  } else if (code.type == 'name') {
    char_bank <- toupper(as.character(aa_info$Name))
  }
  
  return(sample(char_bank, n, replace=T))
}

affinity(strreverse('CASSYSRTGSYEQYF'), 'LLWNGPMAV')
affinity(strreverse('CASSQGLAYEQFF'), 'LLWNGPMAV')
affinity(strreverse('CASSVEGPGELFF'), 'LLWNGPMAV')
affinity('ILGGGCCCQL', 'LLWNGPMAV')

affinity(strreverse('CASSSGQLTNTEAFF'), 'GLCTLVAML')
affinity(strreverse('CASSFGVNSDYTF'), 'KGYVYQGL')

tcr.seq <- 'SADRVGNT'
tcr.seq.pruned <- 'SADRVGNT'
epitope <- 'INFDFNTI'

random_affinities <- c()

for (i in seq(1,1000)) {
  aff <- affinity(
    paste0(random_aa(str_length(tcr.seq.pruned)), collapse=""),
    epitope
  )
  random_affinities <- c(random_affinities, aff$min.rt.energy)
}
hist(random_affinities)
sd(random_affinities)
mean(random_affinities)

actual_aff <- affinity(tcr.seq.pruned, epitope)
actual_aff
abs(mean(random_affinities) - actual_aff$min.rt.energy) / sd(random_affinities)






