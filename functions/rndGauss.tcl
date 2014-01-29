proc rndGauss {mean sd} {
        #
        # generate a random number that is normally distibuted
        # using the central limit theorem
        #
        # mean -> expected mean of the generated numbers
        # sd -> expected standard deviation of the generated numbers
        #
        
        # number of iterations (the higher, the better, the slower):
        set n 150
        # produce n equally random integers from [0,1]
        set sum 0
        for {set i 0} {$i < $n} {incr i} {
                set sum [expr {$sum + int(rand()*2)}]
        }
        # transform the sum to the interval [0,1] again -> sum/n
        # and then transform to [mean,sd^2]
        #
        return [expr {((($sum/double($n))-0.5) * sqrt($n*12)) * $sd + $mean}]
 }