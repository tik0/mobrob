proc normImportanceFactorsLog {} {
global samples_k numSamples
  set alphaLog [getNormFactorLog]
  # set the new importanceFactor
  for {set kk 0} {$kk < $numSamples} {incr kk} {
    # Calculate the real probability
	set importanceFactorLog [ expr [lindex $samples_k($kk) 3] - $alphaLog ]
	if { ($importanceFactorLog == -inf) || ($importanceFactorLog < -500) } {
	   set importanceFactor 0.0
	} else {
	   set importanceFactor [ expr exp($importanceFactorLog) ]
	}
	set samples_k($kk) [list [lindex $samples_k($kk) 0] [lindex $samples_k($kk) 1] [lindex $samples_k($kk) 2] $importanceFactor] 
  }
}