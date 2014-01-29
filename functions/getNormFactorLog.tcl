proc getNormFactorLog {} {
global samples_k numSamples
  set alphaLog -inf
  for {set kk 0} {$kk < $numSamples} {incr kk} {
	set alphaLog [ logAdd $alphaLog [ lindex $samples_k($kk) 3 ] ]
  }
  return $alphaLog
}