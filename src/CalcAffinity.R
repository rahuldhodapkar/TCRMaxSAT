#!/usr/bin/env R
# CalcAffinity.R
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
library(stringr)
library(dplyr)
library(ggplot2)

homo_sapiens_trb_vdjdb <- read.csv('./config/homo_sapiens_trb_vdjdb.tsv', header=T, sep='\t')

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

#' Generate random amino acid sequence of defined length
#'
#' \code{random_cdr3} returns a random amino acid sequence
#'
#' @param n Integer length of the amino acid sequence desired
#' @param code.type string, one of 'olc', 'tlc', or 'name', defining
#'   in what form the amino acid codes will be
#' @return a character vector of an amino acid sequence
#' @examples
#' random_cdr3()
random_cdr3 <- function(n = 7, code.type = 'olc') {
  cdr3_bank <- homo_sapiens_trb_vdjdb$CDR3
  # WARNING "n" not yet supported.
  # WARNING only olc supported currently.
  
  return(as.character(sample(cdr3_bank, 1)))
}

#' Rate interaction of CDR3 against simulated interactions
#'
#' \code{simulate_ixns} runs simulated cdr3-epitope interactions
#'   and ranks provided cdr3 against simulations.
#'
#' @param cdr3 String; AA sequence of the cdr3 region
#' @param epitope String; AA sequence of the epitope
#' @param num.sims Integer; number of simulations to run, default 1000
#' @return mean energy of simulated sequences,
#'   energy difference between tcr sequence and CDR3, as well as
#'   a z-score for this difference
#' @examples
#' simulate_ixns('CASSEATGASYEQYF', 'LLWNGPMAV')
simulate_ixns <- function(cdr3, epitope, num.sims = 1000, aa.generator=random_aa) {
  simulated_rand_affinities <- c()
  
  for (i in seq(1, num.sims)) {
    aff <- affinity(
      paste0(aa.generator(str_length(cdr3)), collapse=""),
      epitope
    )
    simulated_rand_affinities <- c(simulated_rand_affinities, aff$min.rt.energy)
  }
  actual_aff <- affinity(cdr3, epitope)
  
  return(list(
    mean_sim_energy = mean(simulated_rand_affinities),
    cdr3_energy = actual_aff$min.rt.energy,
    sd_sim_energy = sd(simulated_rand_affinities),
    diff_sim_cdr3 = mean(simulated_rand_affinities) - actual_aff$min.rt.energy,
    z.score = (mean(simulated_rand_affinities) - actual_aff$min.rt.energy)/sd(simulated_rand_affinities)
  ))
}

homo_sapiens_trb_vdjdb_hq <-subset(homo_sapiens_trb_vdjdb, Score >= 3)

TRIM_FRONT <- 3;
TRIM_BACK <- 3;
set.seed(103)
mean_sim_energies <- c()
cdr3_energies <- c()
diff_sim_cdr3s <- c()
z_scores <- c()
shuffled_ixs <- sample(seq(1, nrow(homo_sapiens_trb_vdjdb_hq)))

for (i in shuffled_ixs[1:50]) {
  print(i)
  cdr3_seq <- strreverse(as.character(homo_sapiens_trb_vdjdb_hq$CDR3[[i]]))
  cdr3_trimmed <- substr(cdr3_seq, TRIM_FRONT, str_length(cdr3_seq) - TRIM_BACK)
  simout <- simulate_ixns(cdr3_trimmed, 
                          as.character(homo_sapiens_trb_vdjdb_hq$Epitope[[i]]),
                          aa.generator = random_cdr3)
  mean_sim_energies <- c(mean_sim_energies, simout$mean_sim_energy)
  cdr3_energies <- c(cdr3_energies, simout$cdr3_energy)
  diff_sim_cdr3s <- c(diff_sim_cdr3s, simout$diff_sim_cdr3)
  z_scores <- c(z_scores, simout$z.score)
}

hist(mean_sim_energies)
hist(cdr3_energies)
hist(diff_sim_cdr3s)
hist(z_scores)

df <- data.frame(
  energy = c(mean_sim_energies, cdr3_energies),
  source = c(rep("Sim", length(mean_sim_energies)),
            rep("CDR3", length(cdr3_energies)))
)
df$source <- as.factor(df$source)
mu <- df %>% 
        group_by(source) %>%
        summarise(grp.mean=mean(energy))

ggplot(df, aes(x=energy, fill=source)) +
  geom_density(alpha=0.4) + 
  geom_vline(data=mu, aes(xintercept=grp.mean, color=source),
               linetype="dashed") +
  ggtitle("RT Energy of Linearized Binding Model")

ggplot(data.frame(energy = z_scores), aes(x=energy)) + 
  geom_density(alpha=0.4)

ggplot(data.frame(energy = diff_sim_cdr3s), aes(x=energy)) + 
  geom_density(alpha=0.4)

