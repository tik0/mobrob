proc isnan {x} {
    if { ![string is double $x] || $x != $x } {
        return 1
    } else {
        return 0
    }
}
