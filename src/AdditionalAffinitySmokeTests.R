#!/usr/bin/env R
# AdditionalAffinitySmokeTests.R
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

#################################################
## Repeat with decoys
#################################################

TRIM_FRONT <- 3;
TRIM_BACK <- 3;
set.seed(103)
decoy_mean_sim_energies <- c()
decoy_cdr3_energies <- c()
decoy_diff_sim_cdr3s <- c()
decoy_z_scores <- c()
shuffled_ixs <- sample(seq(1, nrow(homo_sapiens_trb_vdjdb_hq)))

for (i in shuffled_ixs[1:50]) {
  print(i)
  cdr3_seq <- paste0(
    random_aa(str_length(as.character(homo_sapiens_trb_vdjdb_hq$CDR3[[i]]))),
    collapse = '')
  cdr3_trimmed <- substr(cdr3_seq, TRIM_FRONT, str_length(cdr3_seq) - TRIM_BACK)
  simout <- simulate_ixns(cdr3_trimmed, 
                          as.character(homo_sapiens_trb_vdjdb_hq$Epitope[[i]]))
  decoy_mean_sim_energies <- c(decoy_mean_sim_energies, simout$mean_sim_energy)
  decoy_cdr3_energies <- c(decoy_cdr3_energies, simout$cdr3_energy)
  decoy_diff_sim_cdr3s <- c(decoy_diff_sim_cdr3s, simout$diff_sim_cdr3)
  decoy_z_scores <- c(decoy_z_scores, simout$z.score)
}

hist(decoy_mean_sim_energies)
hist(decoy_cdr3_energies)
hist(decoy_diff_sim_cdr3s)
hist(decoy_z_scores)

decoy_df <- data.frame(
  energy <- c(decoy_mean_sim_energies, decoy_cdr3_energies),
  source <- c(rep("Sim", length(decoy_mean_sim_energies)),
              rep("Decoy", length(decoy_cdr3_energies)))
)
decoy_df$source <- as.factor(df$source)
decoy_mu <- decoy_df %>% 
  group_by(source) %>%
  summarise(grp.mean=mean(energy))

ggplot(decoy_df, aes(x=energy, fill=source)) +
  geom_density(alpha=0.4) + 
  geom_vline(data=decoy_mu, aes(xintercept=grp.mean, color=source),
             linetype="dashed")

ggplot(data.frame(energy = decoy_z_scores), aes(x=energy)) + 
  geom_density(alpha=0.4)

ggplot(data.frame(energy = decoy_diff_sim_cdr3s), aes(x=energy)) + 
  geom_density(alpha=0.4)

#################################################
## Repeat with non-matched CDR3s
#################################################

TRIM_FRONT <- 3;
TRIM_BACK <- 3;
set.seed(103)
mismatch_mean_sim_energies <- c()
mismatch_cdr3_energies <- c()
mismatch_diff_sim_cdr3s <- c()
mismatch_z_scores <- c()
shuffled_ixs <- sample(seq(1, nrow(homo_sapiens_trb_vdjdb_hq)))

for (i in shuffled_ixs[1:50]) {
  print(i)
  epitope_seq <- homo_sapiens_trb_vdjdb_hq$Epitope[[i]];
  cdr3_seq <- sample(homo_sapiens_trb_vdjdb_hq$CDR3[
    homo_sapiens_trb_vdjdb_hq$Epitope != epitope_seq
    ], 1)
  cdr3_trimmed <- substr(cdr3_seq, TRIM_FRONT, str_length(cdr3_seq) - TRIM_BACK)
  simout <- simulate_ixns(cdr3_trimmed, 
                          as.character(epitope_seq))
  mismatch_mean_sim_energies <- c(mismatch_mean_sim_energies, simout$mean_sim_energy)
  mismatch_cdr3_energies <- c(mismatch_cdr3_energies, simout$cdr3_energy)
  mismatch_diff_sim_cdr3s <- c(mismatch_diff_sim_cdr3s, simout$diff_sim_cdr3)
  mismatch_z_scores <- c(mismatch_z_scores, simout$z.score)
}

hist(mismatch_mean_sim_energies)
hist(mismatch_cdr3_energies)
hist(mismatch_diff_sim_cdr3s)
hist(mismatch_z_scores)

mismatch_df <- data.frame(
  energy <- c(mismatch_mean_sim_energies, mismatch_cdr3_energies),
  source <- c(rep("Sim", length(mismatch_mean_sim_energies)),
              rep("Diff", length(mismatch_cdr3_energies)))
)
mismatch_df$source <- as.factor(mismatch_df$source)
mismatch_mu <- mismatch_df %>% 
  group_by(source) %>%
  summarise(grp.mean=mean(energy))

ggplot(mismatch_df, aes(x=energy, fill=source)) +
  geom_density(alpha=0.4)

ggplot(data.frame(energy = mismatch_z_scores), aes(x=energy)) + 
  geom_density(alpha=0.4)

ggplot(data.frame(energy = mismatch_diff_sim_cdr3s), aes(x=energy)) + 
  geom_density(alpha=0.4)