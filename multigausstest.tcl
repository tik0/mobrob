# ===========================================================================
# 
# multigausstest.tcl --
# test of multivariate Gaussian distribution generator
# 
# Ralf Moeller
# 
#    Copyright (C) 2007
#    Computer Engineering Group
#    Faculty of Technology
#    University of Bielefeld
#    www.ti.uni-bielefeld.de
# 
# 1.0 / 12. Jun 07 (rm)
# - from scratch
# 
# 
# ===========================================================================

#----------------------------------------------------------------------------
# rndMultiVarGauss draws a single sample from a multivariate Gaussian
# distribution specified by a center and a covariance matrix
#----------------------------------------------------------------------------
# C_name is the name of a covariance matrix (input)
# center is the name of a center column vector (input)
# x_name is the name of a randomly drawn column vector (output)
proc rndMultiVarGauss {C_name center_name x_name} {
    # enter new matrix context
    mat enter
    # local variables
    # matrices for eigenvectors, eigenvalues, stddev's, random numbers
    mat matrix E
    mat matrix lambda
    mat matrix sigma
    mat matrix rnd
    # eigenvalue decomposition
    mat eigenSymm $C_name E lambda
    # puts [formatMatrix E]
    # puts [formatMatrix lambda]
    # lambda = sigma^2 -> sigma = sqrt(lambda)
    mat sqrt lambda sigma
    # generate random numbers with stddev. according to eigenvalues
    set rndList {}
    set sigmaList [mat getVectorList sigma]
    foreach s $sigmaList {
	lappend rndList [gsl randGaussian $s]
    }
    mat setRowVector rnd $rndList
    # assemble random vector
    # each row of E is multiplied element-wise by rnd
    mat multRows E rnd E
    # rnd1 * col1 + rnd2 * col2 + ...
    mat sumAlongRows E $x_name
    # add center vector
    mat addTo $center_name $x_name
    # leave matrix context (deletes local variables)
    mat leave
}


# initialize random number generator
gsl randSeed [clock seconds]

# generate C
mat matrix C
# try different examples:
# mat setRowByRow C {{1 0 0} {0 4 0} {0 0 0.5}}
# mat setRowByRow C {\
		       {1.02182 -0.72259 -0.196272} \
		       {-0.72259 0.75134 0.151482} \
		       {-0.196272 0.151482 0.05205}}
mat setRowByRow C {\
		       {5.39595 -7.31075 -0.712898} \
		       {-7.31075 11.4388 0.935503} \
		       {-0.712898 0.935503 0.11295}}

# canvas coordinates and scaling
set leftW1 0.0
set bottomW1 0.0
set widthW1 30.0
set heightW1 30.0
set wToC1 20.0

# define center vector (current robot location)
set x [expr $leftW1 + $widthW1/2]
set y [expr $bottomW1 + $heightW1/2]
set theta 0.0
mat matrix center
mat setColVector center [list $x $y $theta]
# puts "center [mat getVectorList center]"

# create canvas
createWorldCanvas world1 $leftW1 $bottomW1 $widthW1 $heightW1 $wToC1

# draw raudom numbers and visualize them as dots
mat matrix rndVec
for {set t 0} {$t < 1000} {incr t} {
    rndMultiVarGauss C center rndVec
    set rndVecList [mat getVectorList rndVec]
    lassign $rndVecList xr yr thetar
    drawCircleWC world1 $xr $yr 2 darkgrey lightgrey 1 "dot"
    update
}

# draw uncertainty ellipse of C
drawUncertainty world1 $x $y $theta C 1.0 100 2 red "ellipse"

