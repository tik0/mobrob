proc logAdd {a b} {
    set c 0.0
    
    if {$a > $b} {
        set c $b
        set b $a
        set a $c
    }
    set c [expr $a - $b]
    
    if { ($a == -inf) || ($c < -36) } {
        return $b
    } else {
        return [expr $b + log(exp($c) + 1.0 )];
    }
}