proc normImportanceFactorsLog {} {
global samples_k numSamples
  set alphaLog [getNormFactorLog]
  # set the new importanceFactor
  for {set kk 0} {$kk < $numSamples} {incr kk} {
    # Calculate the real probability
    if {[catch {set importanceFactor [expr exp([lindex $samples_k($kk) 3] - $alphaLog)]}]} {
	  set importanceFactor 0.0
	}
	set samples_k($kk) [list [lindex $samples_k($kk) 0] [lindex $samples_k($kk) 1] [lindex $samples_k($kk) 2] $importanceFactor] 
  }
}