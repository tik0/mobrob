proc checkPosEigVal {C_name} {
    # enter new matrix context
    mat enter
    # local variables
    # matrices for eigenvectors, eigenvalues
    mat matrix E
    mat matrix lambda
    # eigenvalue decomposition
    mat eigenSymm $C_name E lambda
    set lambdaList [mat getVectorList lambda]
	mat leave
    foreach s $lambdaList {
		if { $s < 0.0 } {
			return 0
	    }
    }
	return 1
	
}