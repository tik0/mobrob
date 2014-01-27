# ===========================================================================
# 
# A3.tcl --
#
# 0.1
# ===========================================================================

frame .buFrame
# pack .buDraw .buErase .buExit -in .buFrame -side top -anchor n
pack .buFrame -in .
set cWidth 300.0
set cHight 300.0
canvas .c -width $cWidth -height $cHight
pack .c -in .
set particleLenght 10



puts "([info script])"

robosim clearWorld
robosim clearRobot

# what to do on startup
set doObstacleAvoidance 1

# robot with differential steering
# wheel radius
set r 0.5
# wheel distance from center
set l 1.0
# wheel width
set w 0.3
# create left and right wheel
set leftWheelID [robosim addWheel $r $w $l $PIHALF 0.0 1 0]
set rightWheelID [robosim addWheel $r $w $l [expr -$PIHALF] $PI 1 0]
# create body
set l2 [expr $l / 2.0]
set ml2 [expr -$l2]
set ml [expr -$l]
robosim addBodyPart $ml2 $l $l2 $l
robosim addBodyPart $l2 $l $l $l2
robosim addBodyPart $l $l2 $l $ml2
robosim addBodyPart $l $ml2 $l2 $ml
robosim addBodyPart $l2 $ml $ml2 $ml
robosim addBodyPart $ml2 $ml $ml $ml2
robosim addBodyPart $ml $ml2 $ml $l2
robosim addBodyPart $ml $l2 $ml2 $l
robosim addBodyPart 0 0 $l 0

# laser scanner (only for obstacle avoidance)
if {$doObstacleAvoidance} {
    set xR 0.0
    set yR 0.0
    set phiR 0.0
    set sweep [expr $PI]
    set nRays 7
    set range 10.0
    set sigmaSlope 0.0
    set sigmaOffset 0.2
    robosim createLaser $xR $yR $phiR $sweep $nRays $range $sigmaSlope $sigmaOffset
}

# load world
# source "simple.world"

# visualize world
set leftW1 0.0
set bottomW1 0.0
set widthW1 30.0
set heightW1 30.0
set wToC1 20.0
# a key handler is passed to this canvas
createWorldCanvas world1 $leftW1 $bottomW1 $widthW1 $heightW1 $wToC1 \
    keyHandler
setWorldCanvasGrid world1 10

# initialize robot position
set x 10.0
set y 10.0
set theta 0.0

# path integration, see Siegwart/Nourbakhsh p.188
proc updatePose {x y theta vl vr dt} {
    global l
    set sl [expr $dt * $vl]
    set sr [expr $dt * $vr]
    set b [expr 2.0 * $l]
    set slsr2 [expr ($sl + $sr) / 2.0]
    set srsl2b [expr ($sr - $sl) / (2.0 * $b)]
    set x [expr $x + $slsr2 * cos($theta + $srsl2b)]
    set y [expr $y + $slsr2 * sin($theta + $srsl2b)]
    set theta [expr $theta + 2.0 * $srsl2b]
    return "$x $y $theta"
}

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

# motion covariance, see Siegwart/Nourbakhsh p.188
proc motionCovariance {CdeltaName x y theta vl vr dt kl kr} {
    set sl [expr $dt * $vl]
    set sr [expr $dt * $vr]
    mat setRowByRow $CdeltaName \
      "{[expr $kr * abs($sr)] 0.0} \
         {0.0 [expr $kl * abs($sl)]}"
}

#----------------------------------------------------------------------------
# robots control
#----------------------------------------------------------------------------

## Monte Carlo Setup
# Amount of samples
set numSamples 100;
# Initial sample
for {set k 0} {$k < $numSamples} {incr k} {
    # x y theta importanceFactor
    # set samples_k($k) [list [expr rand()*$widthW1] [expr rand()*$heightW1] [expr rand()*2.0*$PI] [expr 1.0/$numSamples]]
     set samples_k($k) [list [expr $x] [expr $y] [expr $theta] [expr 1.0/$numSamples]]
}

# time step
set dt 0.1

# Initial Speed
set vl 0.0
set vr 0.0

set thresholdDist 2.0
set vpos 0.5
set vneg -0.5

# factors for motion covariance
set kl 0.01
set kr 0.01
# initial variance in each direction x,y,theta
set c00 0.1
set c11 0.1
set c22 0.01

# Uncertanty of the laser
set sigmaLaser 10.0

# Laser update every stepSensorUpdate step
set stepSensorUpdate 10

# show uncertainty ellipse
set points 200
set scale 1.0

# frozen uncertainty ellipses
set freeze 0

proc obstacleAvoidanceLoop {} {
    global x y theta dt nRays range thresholdDist kl kr vl vr vpos vneg c00 c11 c22 stopLoop numSamples samples_k widthW1 heightW1 PI particleLenght cWidth cHight points scale freeze sigmaLaser stepSensorUpdate

    ###########################################################
    # Monte-Carlo
    
    # Move the particles
#     set samples_k
#     for {set k 0} {$k < $numSamples} {incr k} {
# 	# x y theta importanceFactor
# 	set samples_k_moved($k) [list [expr rand()*$widthW1] [expr rand()*$heightW1] [expr rand()*2.0*$PI] [expr 1.0/$numSamples]]
#     }
    
    

    
    

    
    set stopLoop 0
    mat enter
    mat matrix H
    mat matrix H1
    
    # Matrix for uncertanties
    mat matrix Cp 3 3
    mat set Cp 0 0 $c00
    mat set Cp 1 1 $c11
    mat set Cp 2 2 $c22
    mat matrix Fp
    mat matrix Fdelta
    mat matrix Cdelta
    # define center vector (current robot location)
#     set x [expr $leftW1 + $widthW1/2]
#     set y [expr $bottomW1 + $heightW1/2]
#     set theta 0.0
    mat matrix center
    mat setColVector center [list 0.0 0.0 0.0]
    mat matrix rndVec
    
    set loopCount 0
    while {!$stopLoop} {
      # draw new robot
      deleteRobot world1
      drawRobot world1 $x $y $theta
      # draw uncertainty ellipse
      # puts [formatMatrix Cp]
      deleteEllipse world1
      drawUncertainty world1 $x $y $theta Cp $scale $points
      if {$freeze} {
          drawUncertainty world1 $x $y $theta Cp $scale $points \
            1 black frozen
          set freeze 0
      }
      update
      # update uncertainty
      motionCovariance Cdelta $x $y $theta $vl $vr $dt $kl $kr
      # puts [formatMatrix Cdelta]
      motionJacobian Fdelta $x $y $theta $vl $vr $dt
      poseJacobian Fp $x $y $theta $vl $vr $dt
      # update Cp from Cp and Cdelta and the Jacobians
      # (error propagation law)
      matset Cp [matexpr Fp * Cp * ~Fp + Fdelta * Cdelta * ~Fdelta]
      # puts [formatMatrix Cp]
      # puts "\{[mat getRowByRow Cp]\}"
      # get transformation matrix for current pose
      
      # get transformation matrix for current pose
      robosim getTransformationMatrix $x $y $theta H
      # get laser scan
      set scan [robosim laserScan H]
      # analyze laser scan
      # find closest distance on the left side
      set leftMinDist $range
      for {set n 0} {$n < $nRays / 2} {incr n} {
          set beam [lindex $scan $n]
          set dist [lindex $beam 4]
          if {$dist < $leftMinDist} {
            set leftMinDist $dist
          }
      }
      # find closest distance on the right side
      set rightMinDist $range
      for {set n [expr $nRays / 2]} {$n < $nRays} {incr n} {
          set beam [lindex $scan $n]
          set dist [lindex $beam 4]
          if {$dist < $rightMinDist} {
            set rightMinDist $dist
          }
      }
      # analyse both sides
      set leftClose [expr $leftMinDist < $thresholdDist]
      set rightClose [expr $rightMinDist < $thresholdDist]
      if {$leftClose && $rightClose} {
          set vl $vneg
          set vr $vneg
      } elseif {$leftClose} {
          set vl $vpos
          set vr $vneg
      } elseif {$rightClose} {
          set vl $vneg
          set vr $vpos
      } else {
          set vl [expr 2 * $vpos]
          set vr [expr 2 * $vpos]
      }
      # update pose
      lassign [updatePose $x $y $theta $vl $vr $dt] x1 y1 theta1
      # get transformation matrix for next pose
      robosim getTransformationMatrix $x1 $y1 $theta1 H1
      # check collision for current position
      set c [robosim collision H]
      if {$c} {
          # there's nothing we can do
          puts "collision at current position"
          break
      }
      # check collision for predicted position and movement
      set c1 [robosim collision H1]
      set mc [robosim movementCollision H H1]
      if {$c1 || $mc} {
          puts "predicted collision at $x1 $y1 $theta1"
          # stop the robot
          set vl 0.0
          set vr 0.0
          # draw a pair of collision dots
          drawCircleWC world1 $x $y 3 orange yellow 1 "collision"
          drawCircleWC world1 $x1 $y1 3 orange yellow 1 "collision"
      } else {
          # no collision: update robot pose
          set x $x1
          set y $y1
          set theta $theta1
      }
      
      
      ## Animation of particle
      # Scaling factors
      set scaleWidth [expr $cWidth / $widthW1]
      set scaleHight [expr $cHight / $heightW1]
      # Draw robot
      .c delete robot
      .c create line [expr $x * $scaleWidth]  [expr $cHight - $y * $scaleHight]  [expr $x*$scaleWidth+$particleLenght*cos($theta)] [expr $cHight - $y*$scaleHight -$particleLenght*sin($theta)] -arrow last -fill black -tags robot
      # Draw particles
      .c delete particles
      for {set k 0} {$k < $numSamples} {incr k} {
		## Importance sampling (draws a sample from samples_k)
		set kk 0
		set randNr [expr rand()]
		set propInt 0 
		for {set kIdx 0} {$kIdx < $numSamples} {incr kIdx} {
			# Draw a sample between 0 .. 1
			# Integrate over importance Factors (index 3)
			set propInt [expr $propInt + [lindex $samples_k($kIdx) 3]]
			# Get kk as the proper index for the sample
			if {$propInt >= $randNr} {
			  set kk $kIdx
			  break
			  }
		  }
		# puts $kk
      
		# Get the particle
		set pxOld [lindex $samples_k($kk) 0]
		set pyOld [lindex $samples_k($kk) 1]
		set pthetaOld [expr [lindex $samples_k($kk) 2]]
		
		# update particle
		lassign [updatePose $pxOld $pyOld $pthetaOld $vl $vr $dt] x1 y1 theta1
		
		## Add noise	
		# Get random values
		mat matrix CpInit 3 3
		#mat set CpInit 0 0 [ expr $c00 * abs($pxOld - $x1)]
		#mat set CpInit 1 1 [ expr $c11 * abs($pyOld - $y1)]
		#mat set CpInit 2 2 [ expr $c22 * abs($pthetaOld - $theta1)]
		
		matset CpInit [matexpr Fdelta * Cdelta * ~Fdelta]
		
	
		
		rndMultiVarGauss CpInit center rndVec
		set rndVecList [mat getVectorList rndVec]
		lassign $rndVecList xr yr thetar
		
		puts "$xr $yr $thetar"
		lassign [isnan $xr ] xrIsNan
		lassign [isnan $yr ] yrIsNan
		lassign [isnan $thetar ] thetarIsNan
		if { $xrIsNan } {
			set xr 0.0
		}
		if { $yrIsNan } {
			set yr 0.0
		}
		if { $thetarIsNan } {
			set thetar 0.0
		}
		
		
		set ptheta [expr $theta1 + $thetar]
		set px [expr $x1 + $xr]
		set py [expr $y1 + $yr]
		
		# Draw
		.c create line [expr $px * $scaleWidth]  [expr $cHight - $py * $scaleHight]  [expr $px*$scaleWidth+$particleLenght*cos($ptheta)] [expr $cHight - $py*$scaleHight -$particleLenght*sin($ptheta)] -arrow last -fill red -tags particles
			

			## get transformation matrix for current particle
			robosim getTransformationMatrix $px $py $ptheta H
			# get laser scan from particle
			set scanParticle [robosim laserScan H]
			set importanceFactor 0.0
			for {set n 0} {$n < $nRays} {incr n} {
				# Get the laser value from the particle
				set beamParticle [lindex $scanParticle $n]
				set distParticle [lindex $beamParticle 4]
				# Get the laser value from the robot
				set beamLaser [lindex $scan $n]
				set distLaser [lindex $beamLaser 4]
				if {$distLaser > $range} {
					set distLaser $range
				}
				if {$distParticle > $range} {
					set distParticle $range
				}
				set importanceFactor [expr $importanceFactor - log($sigmaLaser * sqrt(2.0 * $PI)) - pow(($distLaser - $distParticle),2) / 2.0 / pow($sigmaLaser, 2) ]
			}
			
			set drawList_k($k) $kk
			set importanceFactor_k($k) $importanceFactor
			set px_k($k) $px
			set py_k($k) $py
			set ptheta_k($k) $ptheta
		  }
		  
		  
		incr loopCount
		# Odometryupdate else Laserupdate
		if  { $loopCount % $stepSensorUpdate != 0 } {
			for {set kk 0} {$kk < $numSamples} {incr kk} {
				set samples_k($kk) [list $px_k($kk) $py_k($kk) $ptheta_k($kk) [lindex $samples_k($kk) 3]] 
			}
		} else {

		  ### Assigne new values for particles
		  # Init the new particles
		  for {set kk 0} {$kk < $numSamples} {incr kk} {
	# 	set importanceFactorTmp_k($value) [expr exp($importanceFactor_k($value))]
			set samples_k($kk) [list 0.0 0.0 0.0 0.0] 
		  }
		  
		  # Save the new particles
		  foreach {key value} [array get drawList_k] {
		# 	set importanceFactorTmp_k($value) [expr exp($importanceFactor_k($value))]
			set samples_k($value) [list $px_k($key) $py_k($key) $ptheta_k($key) [expr exp($importanceFactor_k($key))] ] 
			#if {[catch {set samples_k($value) [list $px_k($key) $py_k($key) $ptheta_k($key) [expr exp($importanceFactor_k($key))] ] } errmsg]} {
			#  set samples_k($value) [list $px_k($key) $py_k($key) $ptheta_k($key) -1.0] 
			#}
		  }
		  
		  ## Renormalize 
		  set alpha 0.0
		 
		  for {set kk 0} {$kk < $numSamples} {incr kk} {
			set alpha [ expr $alpha + [ lindex $samples_k($kk) 3 ] ]
		  }
		  # set the new importanceFactor
		  for {set kk 0} {$kk < $numSamples} {incr kk} {
			set kIdx $drawList_k($kk)
			set samples_k($kk) [list [lindex $samples_k($kk) 0] [lindex $samples_k($kk) 1] [lindex $samples_k($kk) 2] [expr [lindex $samples_k($kk) 3] / $alpha]] 
		  }
		  
		  
		  ### Resample (remove samples with 0 importance and duplicate the others by drawing)
		  set numNoImportance 0
		  for {set kk 0} {$kk < $numSamples} {incr kk} {
		  #puts [lindex $samples_k($kk) 3]
			if {[lindex $samples_k($kk) 3] == 0.0} {
			  set samplesNoImportance_k($numNoImportance) $kk
			  incr numNoImportance
			}
		  }
		  puts "numNoImportance $numNoImportance"
		  puts "-----------------------------------------------------------------"
		  
		  # Resample the bad ones (1: copy living samples) (2: random samples) (else: let them die)
		  set resampleMethod 1
		  if {$resampleMethod == 1} {
			  foreach {key value} [array get samplesNoImportance_k] {
				# Draw a sample
				  set kk 0
				  set randNr [expr rand()]
				  set propInt 0.0
				  for {set kIdx 0} {$kIdx < $numSamples} {incr kIdx} {
					# Draw a sample between 0 .. 1
					# Integrate over importance Factors (index 3)
					#puts "-----------------------------------------------------------------"
					#puts $kIdx
					#puts $samples_k($kIdx)
					#puts [lindex $samples_k($kIdx) 3]
					set propInt [expr $propInt + [lindex $samples_k($kIdx) 3]]
					# Get kk as the proper index for the sample
					if {$propInt >= $randNr} {
					  set kk $kIdx
					  break
					}
				  }
				set samples_k($value) [list [lindex $samples_k($kk) 0] [lindex $samples_k($kk) 1] [lindex $samples_k($kk) 2] [lindex $samples_k($kk) 3]] 
				#set samples_k($value) samples_k($kk) 
			  }
		    } elseif {$resampleMethod == 2} {
			  # Get maximum probability
			  set maxProbIdx 0.0
			  for {set kk 0} {$kk < $numSamples} {incr kk} {
			    if {$maxProbIdx < [lindex $samples_k($kk) 3]} {
					set maxProbIdx [lindex $samples_k($kk) 3]
				}
			  }
			  foreach {key value} [array get samplesNoImportance_k] {
				set samples_k($value) [list [expr rand()*$widthW1] [expr rand()*$heightW1] [expr rand()*2.0*$PI] [expr $maxProbIdx/2.0]]
			  }
			}
		  
		  ## Renormalize after resampling
		  set alpha 0.0
		 
		  for {set kk 0} {$kk < $numSamples} {incr kk} {
			set alpha [ expr $alpha + [ lindex $samples_k($kk) 3 ] ]
		  }
		  # set the new importanceFactor
		  for {set kk 0} {$kk < $numSamples} {incr kk} {
			set kIdx $drawList_k($kk)
			set samples_k($kk) [list [lindex $samples_k($kk) 0] [lindex $samples_k($kk) 1] [lindex $samples_k($kk) 2] [expr [lindex $samples_k($kk) 3] / $alpha]] 
		  }
		  
		}
#       sleep 1
    }
    mat leave
    

    
}


#----------------------------------------------------------------------------
# rndMultiVarGauss draws a single sample from a multivariate Gaussian
# distribution specified by a center and a covariance matrix
#----------------------------------------------------------------------------
# C_name is the name of a covariance matrix (input)
# center is the name of a center column vector (input)
# x_name is the name of a randomly drawn column vector (output)
proc rndMultiVarGauss {C_name center_name x_name} {
    # enter new matrix context
    mat enter
    # local variables
    # matrices for eigenvectors, eigenvalues, stddev's, random numbers
    mat matrix E
    mat matrix lambda
    mat matrix sigma
    mat matrix rnd
    # eigenvalue decomposition
    mat eigenSymm $C_name E lambda
    # puts [formatMatrix E]
    # puts [formatMatrix lambda]
    # lambda = sigma^2 -> sigma = sqrt(lambda)
    mat sqrt lambda sigma
    # generate random numbers with stddev. according to eigenvalues
    set rndList {}
    set sigmaList [mat getVectorList sigma]
    foreach s $sigmaList {
	lappend rndList [gsl randGaussian $s]
    }
    mat setRowVector rnd $rndList
    # assemble random vector
    # each row of E is multiplied element-wise by rnd
    mat multRows E rnd E
    # rnd1 * col1 + rnd2 * col2 + ...
    mat sumAlongRows E $x_name
    # add center vector
    mat addTo $center_name $x_name
    # leave matrix context (deletes local variables)
    mat leave
}

proc isnan {x} {
    if { ![string is double $x] || $x != $x } {
        return 1
    } else {
        return 0
    }
}

# uncomment these lines for obstacle avoidance movement on startup:
if {$doObstacleAvoidance} {
   #  source "obstacleAvoidance.world"
    source "TestWorld.world"
    drawWorld world1
    obstacleAvoidanceLoop
    
}
