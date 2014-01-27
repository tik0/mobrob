# path integration, see Siegwart/Nourbakhsh p.188
proc updatePose {x y theta vl vr dt} {
    global l
    set sl [expr $dt * $vl]
    set sr [expr $dt * $vr]
    set b [expr 2.0 * $l]
    set slsr2 [expr ($sl + $sr) / 2.0]
    set srsl2b [expr ($sr - $sl) / (2.0 * $b)]
    set x [expr $x + $slsr2 * cos($theta + $srsl2b)]
    set y [expr $y + $slsr2 * sin($theta + $srsl2b)]
    set theta [expr $theta + 2.0 * $srsl2b]
    return "$x $y $theta"
}
