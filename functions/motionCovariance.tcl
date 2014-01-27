# motion covariance, see Siegwart/Nourbakhsh p.188
proc motionCovariance {CdeltaName x y theta vl vr dt kl kr} {
    set sl [expr $dt * $vl]
    set sr [expr $dt * $vr]
    mat setRowByRow $CdeltaName \
      "{[expr $kr * abs($sr)] 0.0} \
         {0.0 [expr $kl * abs($sl)]}"
}
