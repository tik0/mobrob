# ===========================================================================
# 
# A3.tcl --
#
# 0.1
# ===========================================================================

# ### Create a main window for the particles 
frame .buFrame
# pack .buDraw .buErase .buExit -in .buFrame -side top -anchor n
pack .buFrame -in .
set cWidth 300.0
set cHight 300.0
canvas .c -width $cWidth -height $cHight
pack .c -in .
set particleLenght 10

# ### Get some functions
source "functions/rndMultiVarGauss.tcl"
source "functions/isnan.tcl"
source "functions/motionCovariance.tcl"
source "functions/motionJacobian.tcl"
source "functions/updatePose.tcl"
source "functions/drawSample.tcl"
source "functions/poseJacobian.tcl"
source "functions/getNormFactor.tcl"
source "functions/normImportanceFactors.tcl"
source "functions/logAdd.tcl"
source "functions/normImportanceFactorsLog.tcl"
source "functions/getNormFactorLog.tcl"


# ### Common script
puts "([info script])"

robosim clearWorld
robosim clearRobot

# ## robot with differential steering
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

# Uncertanty of the laser
set sigmaLaser 0.5

# ## laser scanner
set xR 0.0
set yR 0.0
set phiR 0.0
set sweep [expr $PI]
set nRays 21
set range 10.0
set sigmaSlope 0.0
set sigmaOffset $sigmaLaser
robosim createLaser $xR $yR $phiR $sweep $nRays $range $sigmaSlope $sigmaOffset

# ## visualize world
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
set x 5.0
set y 12.0
set theta 0.0

# Scaling factors for animation
set scaleWidth [expr $cWidth / $widthW1]
set scaleHight [expr $cHight / $heightW1]


#----------------------------------------------------------------------------
# robots control
#----------------------------------------------------------------------------

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

# Laser update every stepSensorUpdate step
set stepSensorUpdate 100

# show uncertainty ellipse
set points 200
set scale 1.0

# frozen uncertainty ellipses
set freeze 0


## Monte Carlo Setup
# Amount of samples
set numSamples 300;
# Initial sample
for {set k 0} {$k < $numSamples} {incr k} {
    # x y theta importanceFactor
	# Random
    # set samples_k($k) [list [expr rand()*$widthW1] [expr rand()*$heightW1] [expr rand()*2.0*$PI] [expr 1.0/$numSamples]]
	# At robots position
    # set samples_k($k) [list [expr $x] [expr $y] [expr $theta] [expr 1.0/$numSamples]]
	# At robots position with noise
	set samples_k($k) [list [expr $x+ [gsl randGaussian $c00]] [expr $y + [gsl randGaussian $c11]] [expr $theta + [gsl randGaussian $c22]] [expr 1.0/$numSamples]]
}

proc obstacleAvoidanceLoop {} {
    global x y theta dt nRays range thresholdDist kl kr vl vr vpos vneg c00 c11 c22 stopLoop numSamples samples_k widthW1 heightW1 PI particleLenght cWidth cHight points scale freeze sigmaLaser stepSensorUpdate scaleWidth scaleHight
    
    set stopLoop 0
    mat enter
    mat matrix H
    mat matrix H1
    
    # Matrix for uncertanties
    mat matrix Cp 3 3
    mat set Cp 0 0 $c00
    mat set Cp 1 1 $c11
    mat set Cp 2 2 $c22
    matset CpInit [matexpr Cp]
	
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
      #deleteEllipse world1
      #drawUncertainty world1 $x $y $theta Cp $scale $points
      #if {$freeze} {
      #    drawUncertainty world1 $x $y $theta Cp $scale $points \
      #      1 black frozen
      #    set freeze 0
      #}
      #update
      # update uncertainty
      motionCovariance Cdelta $x $y $theta $vl $vr $dt $kl $kr
      # puts [formatMatrix Cdelta]
      motionJacobian Fdelta $x $y $theta $vl $vr $dt
      poseJacobian Fp $x $y $theta $vl $vr $dt
      # update Cp from Cp and Cdelta and the Jacobians
      # (error propagation law)
      matset Cp [matexpr Fp * Cp * ~Fp + Fdelta * Cdelta * ~Fdelta]
      matset Tmp [matexpr Fdelta * Cdelta * ~Fdelta]
      # puts [formatMatrix Cp]
      # puts "\{[mat getRowByRow Cp]\}"
      # get transformation matrix for current pose
      
      # get transformation matrix for current pose
      robosim getTransformationMatrix $x $y $theta H
      # get laser scan
      set scan [robosim laserScan H]
      drawLaserScan world1 $scan
      update
      # analyse laser scan
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
	  set vlOld $vl
      set vrOld $vr
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
      
      
      # ### Calculation of the particle movement

      # Draw particles
      .c delete particles
      for {set k 0} {$k < $numSamples} {incr k} {
		## Importance sampling only when we have a sensor update (draws a sample from samples_k)
		if  { $loopCount % $stepSensorUpdate != 0 } {
			set kk $k
		} else {
			set kk [drawSample]
		}
      
		# Get the old particle
		set pxOld [lindex $samples_k($kk) 0]
		set pyOld [lindex $samples_k($kk) 1]
		set pthetaOld [expr [lindex $samples_k($kk) 2]]
		
		# Update particle with odometry
		lassign [updatePose $pxOld $pyOld $pthetaOld $vl $vr $dt] x1 y1 theta1

		if  { $loopCount % $stepSensorUpdate != 0 } {
			set ptheta $theta1
			set px $x1
			set py $y1
		} else {
			# ## Add noise from Gaussian distribution
			# Get random values from steering uncertainties
			# update uncertainty
			motionCovariance Cdelta $pxOld $pyOld $pthetaOld $vlOld $vrOld $dt $kl $kr
			# puts [formatMatrix Cdelta]
			motionJacobian Fdelta $pxOld $pyOld $pthetaOld $vlOld $vrOld $dt
			poseJacobian Fp $pxOld $pyOld $pthetaOld $vlOld $vrOld $dt

			mat matrix pCp 3 3
			# matset pCp [matexpr Fdelta * Cdelta * ~Fdelta]
			matset pCp [matexpr Fp * CpInit * ~Fp + Fdelta * Cdelta * ~Fdelta]
			rndMultiVarGauss pCp center rndVec
			set rndVecList [mat getVectorList rndVec]
			lassign $rndVecList xr yr thetar
					
			# Add the noise
			set ptheta [expr $theta1 + $thetar]
			set px [expr $x1 + $xr]
			set py [expr $y1 + $yr]
		}
		
		# ## Animation of particles
		.c create line [expr $px * $scaleWidth]  [expr $cHight - $py * $scaleHight]  [expr $px*$scaleWidth+$particleLenght*cos($ptheta)] [expr $cHight - $py*$scaleHight -$particleLenght*sin($ptheta)] -arrow last -fill red -tags particles
			
		# ## Save the new particles
		set px_k($k) $px
		set py_k($k) $py
		set ptheta_k($k) $ptheta
	  }
	   # Draw robot
      .c delete robot
      .c create line [expr $x * $scaleWidth]  [expr $cHight - $y * $scaleHight]  [expr $x*$scaleWidth+$particleLenght*cos($theta)] [expr $cHight - $y*$scaleHight -$particleLenght*sin($theta)] -arrow last -fill black -tags robot
	  
	# Draw the hypotheses of the robot
	  deleteEllipse world1
	  # Calculate the Gaussian out of the samples
	  # get mean (Circular mean is wrong, there we need something imaginary)
	  set mu_x 0.0
	  set mu_y 0.0
	  set mu_t 0.0
	  for {set k 0} {$k < $numSamples} {incr k} {
			set mu_x [expr $mu_x + [lindex $px_k($k)]]
			set mu_y [expr $mu_y + [lindex $py_k($k)]]
			set mu_t [expr $mu_t + [lindex $ptheta_k($k)]]
		}
	  set mu_x [expr $mu_x / $numSamples]
      set mu_y [expr $mu_y / $numSamples]
	  set mu_t [expr $mu_t / $numSamples]
	  # get covariance
	  mat matrix X 3 $numSamples
	  mat matrix Norm 3 3
	  mat fill Norm [expr 1.0 / $numSamples]
	  for {set k 0} {$k < $numSamples} {incr k} {
		  mat set X 0 $k [expr $px_k($k) - $mu_x]
		  mat set X 1 $k [expr $py_k($k) - $mu_y]
		  mat set X 2 $k [expr $ptheta_k($k) - $mu_t]
		}
		mat matrix Sigma
		matset Sigma [matexpr X * ~X]
		mat multE Sigma Norm Sigma
      drawUncertainty world1 $mu_x $mu_y $mu_t Sigma $scale $points
	  update		  
		  
		
		# Odometryupdate else Laserupdate
		incr loopCount
		if  { $loopCount % $stepSensorUpdate != 0 } {
		# #####################################################################
		# #                   Odometry                                        #
		# #####################################################################
			for {set kk 0} {$kk < $numSamples} {incr kk} {
				set samples_k($kk) [list $px_k($kk) $py_k($kk) $ptheta_k($kk) [lindex $samples_k($kk) 3]] 
			}
			puts "Step: $loopCount -> Odometry"
		} else {
		# #####################################################################
		# #                   Laserupdate                                     #
		# #####################################################################
		puts "Step: $loopCount -> Laserupdate"
		
		# ### Get the importance factors
		  for {set kk 0} {$kk < $numSamples} {incr kk} {
			# ## get transformation matrix for current particle
			robosim getTransformationMatrix $px_k($kk) $py_k($kk) $ptheta_k($kk) H
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
				# Calculate the LOGARITHMIC probability log(MUL_n(s_n,l_kk))
				set importanceFactor [expr $importanceFactor - log($sigmaLaser * sqrt(2.0 * $PI)) - pow(($distLaser - $distParticle),2) / 2.0 / pow($sigmaLaser, 2) ]
			}
			# Get P(s,l_kk)
			
			set importanceFactorLog_k($kk) $importanceFactor
			
			#if {[catch {set importanceFactor_k($kk) [expr exp($importanceFactor)]}]} {
			#  set importanceFactor_k($kk) 0.0
			#}
		  }
		  
		  
		  
		  # ### Assign new values for particles
		  # Init the new particles
		  for {set kk 0} {$kk < $numSamples} {incr kk} {
			# set samples_k($kk) [list $px_k($kk) $py_k($kk) $ptheta_k($kk) $importanceFactor_k($kk)]
			set samples_k($kk) [list $px_k($kk) $py_k($kk) $ptheta_k($kk) $importanceFactorLog_k($kk)]
		  }
		  
		  # ## Renormalize samples_k
		  # normImportanceFactors
		  normImportanceFactorsLog
		  # Normalize again (just to be sure in the case, that the log normalization does not work properly)
		  normImportanceFactors
		  }
		  # Check the particle mean
		  set xMean 0.0
		  set yMean 0.0
		  set thetaMean 0.0
		  for {set kk 0} {$kk < $numSamples} {incr kk} {
		   set xMean [expr $xMean + [lindex $samples_k($kk) 0]]
		   set yMean [expr $yMean + [lindex $samples_k($kk) 1]]
		   set thetaMean [expr $thetaMean + [lindex $samples_k($kk) 2]]
		  }
		  set xMean [ expr $xMean / $numSamples]
		  set yMean [ expr $yMean / $numSamples]
		  set thetaMean [ expr $thetaMean / $numSamples]
		  
		  puts "x-Distance [ expr abs($xMean - $x) ]"
		  puts "y-Distance [ expr abs($yMean - $y) ]"
		  puts "theta-Distance [ expr abs($thetaMean - $theta) ]"
		  # sleep 1

		
    }
    mat leave
    

    
}


# Do obstacle avoidance on startup:
# source "obstacleAvoidance.world"
source "TestWorld.world"
drawWorld world1
obstacleAvoidanceLoop
