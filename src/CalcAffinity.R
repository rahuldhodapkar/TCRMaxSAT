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

#' Affinity score between amino acid sequenes.
#'
#' \code{affinity} returns an affinity score between affinity acid sequences
#'
#' @param seq.1 character vector; first amino acid sequence
#' @param seq.2 character vector; second amino acid sequence
#' @param code.type string, one of 'olc', 'tlc', or 'name', defining
#'   in what form the amino acid codes will be in \code{sequence}
#' @return an affinity score that is determined by the residue-residue
#'   energies of amino acids in the sequence
#' @examples
#' affinity(c('I', 'L', 'G', 'G'), c('G', 'G', 'G', 'G'))
#'
#' \dontrun{
#' affinity('ILGG', 'GGGG')
#' }
affinity <- function(seq.1, seq.2, code.type='olc') {
  return(0);
}




