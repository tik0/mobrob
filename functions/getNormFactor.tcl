proc getNormFactor {} {
global samples_k numSamples
  set alpha 0.0
  for {set kk 0} {$kk < $numSamples} {incr kk} {
	set alpha [ expr $alpha + [ lindex $samples_k($kk) 3 ] ]
  }
  return $alpha
}