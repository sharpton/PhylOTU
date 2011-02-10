#tree_to_matrix.R - convert a phylogeny to an R distance matrix
#Copyright (C) 2011  Thomas J. Sharpton 
#author contact: thomas.sharpton@gladstone.ucsf.edu
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#    
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#    
#You should have received a copy of the GNU General Public License
#along with this program (see LICENSE.txt).  If not, see 
#<http://www.gnu.org/licenses/>.

# Invoke % R --slave --args intree outmatrix < tree_to_matrix.R

library(ape)
Args      <- commandArgs()
intree    <- Args[4]
outmatrix <- Args[5]
#refnames  <- Args[6]
phy       <- read.tree( intree )
#reftips   <- as.vector(scan(file=refnames))
#pruned    <- drop.tip(phy, reftips)
#phy$tip.label[grep(pattern, invert=TRUE)])
#phydist   <- cophenetic( pruned )
phydist   <- cophenetic( phy )
write.table( phydist, outmatrix )
