proc drawSample {} {
global samples_k numSamples
 # Draw a sample
 set randNr [expr rand()]
 set propInt 0.0
 for {set kIdx 0} {$kIdx < $numSamples} {incr kIdx} {
  #puts $kIdx
  set propInt [expr $propInt + [lindex $samples_k($kIdx) 3]]
  # Get kk as the proper index for the sample
  if {$propInt >= $randNr} {
   return $kIdx
  }
 }
 # return the last index
 return [expr $numSamples - 1]
}