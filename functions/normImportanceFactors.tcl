proc normImportanceFactors {} {
global samples_k numSamples
  set alpha [getNormFactor]
  # set the new importanceFactor
  for {set kk 0} {$kk < $numSamples} {incr kk} {
	set samples_k($kk) [list [lindex $samples_k($kk) 0] [lindex $samples_k($kk) 1] [lindex $samples_k($kk) 2] [expr [lindex $samples_k($kk) 3] / $alpha]] 
  }
}