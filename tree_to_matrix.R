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
