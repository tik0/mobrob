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
