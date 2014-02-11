# Normalize the angle to [-Pi , PI]
proc sgn {a} {expr {$a>0 ? 1 : $a<0 ? -1 : 0}} ;# rmax

proc normAngle {angle} {
	set PI 3.14159265359
    # Get the whole 2*PI
    set angle_Whole_2PI [ expr floor(abs($angle / 2.0 / $PI)) ]
    # puts $angle_Whole_2PI
    set angle_norm [ expr $angle - $angle_Whole_2PI * 2.0 * $PI]
	if { [expr abs($angle_norm)] > [expr $PI] } {
		set angle_norm [ expr $angle_norm -  [ sgn $angle_norm ] * 2.0 * $PI]
	}
	return $angle_norm
	
}