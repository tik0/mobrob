# ===========================================================================
# 
# A4.tcl --
#
# 0.1
# ===========================================================================

# ### Set variables
set sigmaP 1
set sigmaK 1
set numBecons 10

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
source "functions/rndGauss.tcl"

# ### Common script
puts "([info script])"

robosim clearWorld
robosim clearRobot

# what to do on startup
set doManualControl 1
set doObstacleAvoidance 0

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
    set sweep $PI
    set nRays 21
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

# Initiate beacons
for {set k 0} {$k < $numBecons} {incr k} {
    # x y
	# Random
    set beacon_k($k) [list [expr rand()*$widthW1] [expr rand()*$heightW1]]
	# draw a pair of becons
    drawCircleWC world1 [lindex $beacon_k($k) 0] [lindex $beacon_k($k) 1] 3 red yellow 1 "becons"
}

#----------------------------------------------------------------------------
# robot is controlled manually
#----------------------------------------------------------------------------

# time step
set dt 0.1
# speed state
set vl 0.0
set vr 0.0
# rotatory/translatory speed increase/decrease
set dvr 0.05
set dvt 0.1
# used for keyboard control
set stopLoop 0
# factors for motion covariance
set kl 0.01
set kr 0.01
# initial variance in each direction
set c00 0.1
set c11 0.1
set c22 0.01

proc keyHandler {canvasName key} {
    global u vl vr dvr dvt stopLoop freeze

    switch $key {
      "Left" {
          set vl [expr $vl - $dvr]
          set vr [expr $vr + $dvr]
      }
      "Right" {
          set vl [expr $vl + $dvr]
          set vr [expr $vr - $dvr]
      }
      "Up" {
          set vl [expr $vl + $dvt]
          set vr [expr $vr + $dvt]
      }
      "Down" {
          set vl [expr $vl - $dvt]
          set vr [expr $vr - $dvt]
      }
      "space" {
          set vl 0.0
          set vr 0.0
      }
      "Escape" {
          set stopLoop 1
      }
      "Tab" {
          set freeze 1
      }
    }
}

# show uncertainty ellipse
set points 200
set scale 1.0

# frozen uncertainty ellipses
set freeze 0

set beta 0.0

proc manualControlLoop {} {
    global x y theta dt vl vr stopLoop kl kr c00 c11 c22 scale points freeze beacon_k numBecons sigmaK sigmaP PI beta

	
	
    # source simple.world
    set stopLoop 0
    mat enter
    mat matrix H
    mat matrix H1
    mat matrix Cp 3 3
    mat set Cp 0 0 $c00
    mat set Cp 1 1 $c11
    mat set Cp 2 2 $c22
    mat matrix Fp
    mat matrix Fdelta
    mat matrix Cdelta
    mat matrix pFp
    mat matrix pFdelta
    mat matrix pCdelta
    
    
    # Set the inital covariance and position
    matset Cp_k_k [matexpr Cp]
    
    #matset p_k1_k 3 1
    #matset p_k_k 3 1
    #mat matrix p_k_k
    #mat matrix p_k1_k
    #set pkkList [list $x $y $theta]
    #mat setRowVector pkk $pkkList
    mat matrix p_k1_k 3 1
    mat matrix p_k_k 3 1
    mat set p_k_k 0 0 $x
    mat set p_k_k 1 0 $y
    mat set p_k_k 2 0 $theta
    
    
    while {!$stopLoop} {
    # Set the values
    set p_k_k_List [mat getVectorList p_k_k]

    set px [lindex $p_k_k_List 0]
    set py [lindex $p_k_k_List 1]
    set ptheta [lindex $p_k_k_List 2]
	

	
      
	
      # draw new robot
      deleteRobot world1
      drawRobot world1 $x $y $theta
      # draw uncertainty ellipse
      # puts [formatMatrix Cp]
      deleteEllipse world1
      drawUncertainty world1 $px $py $ptheta Cp_k_k $scale $points
      if {$freeze} {
          drawUncertainty world1 $px $py $ptheta Cp_k_k $scale $points \
            1 black frozen
          set freeze 0
      }
      update
      
      #############################################################################
      # 1. Predict the movement (calculation of p_k1_k and Cp_k1_k
      #############################################################################
      # update uncertainty
      motionCovariance Cdelta $x $y $theta $vl $vr $dt $kl $kr
      # puts [formatMatrix Cdelta]
      motionJacobian Fdelta $x $y $theta $vl $vr $dt
      poseJacobian Fp $x $y $theta $vl $vr $dt
      # update Cp from Cp and Cdelta and the Jacobians
      # (error propagation law)
      
      
      motionCovariance pCdelta $px $py $ptheta $vl $vr $dt $kl $kr
      motionJacobian pFdelta $px $py $ptheta $vl $vr $dt
      poseJacobian pFp $px $py $ptheta $vl $vr $dt
      matset Cp_k1_k [matexpr pFp * Cp_k_k * ~pFp + pFdelta * pCdelta * ~pFdelta]
       # get transformation matrix for current pose
      robosim getTransformationMatrix $px $py $ptheta H
      # update pose
      lassign [updatePose $px $py $ptheta $vl $vr $dt] px1 py1 ptheta1
      # get transformation matrix for next pose
      robosim getTransformationMatrix $px1 $py1 $ptheta1 H1
      
      
      # puts [formatMatrix Cp]
      # puts "\{[mat getRowByRow Cp]\}"
      # get transformation matrix for current pose
      robosim getTransformationMatrix $x $y $theta H
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
          puts "predicted collision at $x $y $theta $vl $vr"
          # stop the robot
          set vl 0.0
          set vr 0.0
          # draw a pair of collision dots
          drawCircleWC world1 $x $y 3 orange yellow 1 "collision"
          drawCircleWC world1 $x1 $y1 3 orange yellow 1 "collision"
      } else {
          # draw line
          drawLineWC world1 $x $y $x1 $y1 blue 1 "trajectory"
          # no collision: update robot pose
          set x $x1
          set y $y1
          set theta $theta1
      }
	# Normalize theta
	#puts $theta
	#if { $theta > 0.0} {
	#  set thetaWhole [ expr floor(abs($theta / 2.0 / $PI)) ]
	#} else {
	#  set thetaWhole [ expr floor(abs($theta / 2.0 / $PI)) ]
	#}
	## puts $thetaWhole
	#set theta [ expr $theta - $thetaWhole * 2.0 * $PI - $PI]
	#puts "Theta $theta"

      
    mat set p_k1_k 0 0 $px1
    mat set p_k1_k 1 0 $py1
    mat set p_k1_k 2 0 $ptheta1
      
      #############################################################################
      # 2. Read sensordata from one arb. beacon (calculation of alpha=z_tilde_P_k1 and beta=z_tilde_K_k1)
      #############################################################################
      # [lindex $p_k_k 0] [lindex $p_k_k 1] [lindex $p_k_k 2] [lindex $u 0] [lindex $u 1] 
      # Get random becon
	set K [ expr int( rand() * $numBecons ) ]
	set atanTmp [expr atan2([lindex $beacon_k($K) 1] - $y, [lindex $beacon_k($K) 0] - $x)]
	set alpha [ expr $atanTmp - $theta - [ rndGauss 0.0 $sigmaP ] ]
	#puts "beacon: $beacon_k($K)"
	#puts "becon: $atanTmp"
	#puts "compas: $thetaCompas"
	#puts " alpha_K $alpha_K"
	#puts $alpha_K
	
	# Get compas
	set beta [ expr atan2($y, $x) + [ rndGauss 0.0 $sigmaK ] ]
	# set alphaNoise [ expr $theta + [ rndGauss 0.0 $sigmaP ] ]
	# puts "Compas: $beta"
	# puts "Compas + Noise: $betaNoise"
	
      #############################################################################
      # 3. Calculate the hypotheses for every beacon i (Calculate z_hat_k1)
      #############################################################################
      
      # Iterate through every beacon
      set matches 0
      for {set i 0} {$i < $numBecons} {incr i} {
          # Get the hypothesis
          set z_hat_k1 [ expr atan2([lindex $beacon_k($i) 1] - $y, [lindex $beacon_k($i) 0] - $x) - $beta]
          
          #############################################################################
	  # 4. Calculate innovation
	  #############################################################################
	  # Get Jacobi-Matrix H for the messurement i at current position
	  set bx [expr [lindex $beacon_k($i) 0] - $x]
	  set by [expr [lindex $beacon_k($i) 1] - $y]
	  set tmp [expr 1.0 / (1.0 + pow($by / $bx, 2))]
	  
	  set h1 [expr $by / pow($bx,2) * $tmp]
	  set h2 [expr -1.0 / $bx * $tmp]
	  set h3 [expr -1.0]
	  mat matrix Hi 1 3
	  mat set Hi 0 0 $h1
	  mat set Hi 0 1 $h2
	  mat set Hi 0 2 $h3
	  
	  # Get covariance R of messurement i
	  set R [expr pow($sigmaP,2)]
	  
	  # Get covariance Cv of innovation
	  matset Cv_tmp [matexpr Hi * Cp_k1_k * ~Hi]
	  
	  set Cv_tmp_list [mat getVectorList Cv_tmp]
	  
	  
	  set Cv [expr [lindex $Cv_tmp_list 0] + $R]
	  
	  # Get innovation v
	  set v [expr $beta - $z_hat_k1]
	  
	  # Calculate the innovation gate and the Kalman gain K
	  set g 1
	  # puts "[expr $v * $Cv * $v]"
	  if { [expr $v * $Cv * $v] <= [expr pow($g,2)] } {
            
	    if {$matches == 0} {
	      incr matches
	      # create matrices for v and Cv
	      mat matrix Cv_mat 1 1
	      mat set Cv_mat 0 0 [expr 1.0 / $Cv]
	      mat matrix v_mat 1 1
	      mat set v_mat 0 0 $v
	  
	      matset K [matexpr Cp_k1_k * ~Hi * Cv_mat]
	      # puts [formatMatrix Hi]
	      
	    } else {
	      incr matches
	      break
	    }
	  }
      }
      #############################################################################
      # 5. Calculate the new position and variance
      #############################################################################
	  
	  #puts "$matches"
      if { $matches == 1} {
         puts "Match found"
         #puts [formatMatrix K]
         matset p_k_k [matexpr p_k1_k + K * v_mat]
         matset Cp_k_k [matexpr Cp_k1_k - K * Cv_mat * ~K]
         
      } else {
	#puts "To many/No matches found"
         matset Cp_k_k [matexpr Cp_k1_k]
         matset p_k_k [matexpr p_k1_k]
      }
     
    }
    mat leave
}

# uncomment this line for manual control on startup:
if {$doManualControl} {
    manualControlLoop
}

#----------------------------------------------------------------------------
# robot is doing a simple obstacle avoidance
#----------------------------------------------------------------------------

set thresholdDist 2.0
set vpos 0.5
set vneg -0.5

proc obstacleAvoidanceLoop {} {
    global x y theta dt nRays range thresholdDist vpos vneg stopLoop
	


    set stopLoop 0
    mat enter
    mat matrix H
    mat matrix H1
    while {!$stopLoop} {
      # draw new robot
      deleteRobot world1
      drawRobot world1 $x $y $theta
      update
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
    }
    mat leave
}

# uncomment these lines for obstacle avoidance movement on startup:
if {$doObstacleAvoidance} {
    source "obstacleAvoidance.world"
    drawWorld world1
    obstacleAvoidanceLoop
}
