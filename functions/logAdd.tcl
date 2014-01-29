proc logAdd {a b} {
    set c 0.0
    
    if {$a > $b} {
        set c $b
        set b $a
        set a $c
    }
    
    if { ($a == -inf) || ([ expr $a - $b ] < -36) } {
        return $b
    } else {
	    set c [expr $a - $b]
        return [expr $b + log(exp($c) + 1.0 )];
    }
}