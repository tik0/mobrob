# pose Jacobian, see Siegwart/Nourbakhsh p.189
proc poseJacobian {FpName x y theta vl vr dt} {
    global l
    set sl [expr $dt * $vl]
    set sr [expr $dt * $vr]
    set b [expr 2.0 * $l]
    set ds [expr ($sl + $sr) / 2.0]
    set dtheta [expr ($sr - $sl) / $b]
    set tt [expr $theta + $dtheta / 2.0]
    set ss [expr $ds / $b]
    set c [expr cos($tt)]
    set s [expr sin($tt)]
    mat setRowByRow $FpName \
      "{1.0 0.0 [expr -$ds * $s]} \
         {0.0 1.0 [expr  $ds * $c]} \
         {0.0 0.0 1.0}"
}
