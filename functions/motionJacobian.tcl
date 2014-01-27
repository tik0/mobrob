# motion Jacobian, see Siegwart/Nourbakhsh p.189
proc motionJacobian {FdeltaName x y theta vl vr dt} {
    global l
    set sl [expr $dt * $vl]
    set sr [expr $dt * $vr]
    set b [expr 2.0 * $l]
    set ds [expr ($sl + $sr) / 2.0]
    set dtheta [expr ($sr - $sl) / $b]
    set tt [expr $theta + $dtheta / 2.0]
    set ss [expr $ds / $b]
    set c [expr 0.5 * cos($tt)]
    set s [expr 0.5 * sin($tt)]
    mat setRowByRow $FdeltaName \
      "{[expr $c - $ss * $s] [expr $c + $ss * $s]} \
         {[expr $s + $ss * $c] [expr $s - $ss * $c]} \
         {[expr 1.0 / $b]      [expr -1.0 / $b]}"
}
