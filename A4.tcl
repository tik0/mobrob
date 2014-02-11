# ===========================================================================
#
# A4.tcl --
#
# 0.1
# ===========================================================================

# ### Set variables
# Standard div. of beacon tracking in rad
set sigmaP 0.01
# Standard div. of compass tracking in rad
set sigmaK 0.01
set numBeacons 2

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
source "functions/checkPosEigVal.tcl"
source "functions/normAngle.tcl"


# ### Common script
# puts "([info script])"

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
for {set k 0} {$k < $numBeacons} {incr k} {
    # x y
    # Random
    set beacon_k($k) [list [expr $widthW1/2.0 + sin(2.0 * $PI * $k / $numBeacons)*$widthW1/2.0] [expr $heightW1/2.0 + cos(2.0 * $PI * $k / $numBeacons)*$heightW1 / 2.0]]
    # set beacon_k($k) [list [expr rand()*$widthW1] [expr rand()*$heightW1]]
	# set beacon_k($k) [list [expr $widthW1/2.0] [expr $heightW1/2.0]]
    # draw a pair of beacons
    drawCircleWC world1 [lindex $beacon_k($k) 0] [lindex $beacon_k($k) 1] 3 red yellow 1 "beacons"
}

#----------------------------------------------------------------------------
# robot is controlled manually
#----------------------------------------------------------------------------

# time step
set dt 0.1
# speed state
set vl 0.5
set vr 0.25
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
    global x y theta dt vl vr stopLoop kl kr c00 c11 c22 scale points freeze beacon_k numBeacons sigmaK sigmaP PI beta



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
	
	
	deleteTag world1 "mean"
	drawCircleWC world1 $px $py 3 orange green 1 "mean"

    # OUTPUT
	set theta [normAngle $theta]
	set ptheta [normAngle $ptheta]
    # puts "theta $theta"
    # puts "ptheta $ptheta"



      # draw new robot
      deleteRobot world1
      drawRobot world1 $x $y $theta
      # draw uncertainty ellipse
      # # puts [formatMatrix Cp]
      deleteEllipse world1
      # puts "-----------------------------"
      # puts "$px $py $ptheta $scale $points"
      # puts [formatMatrix Cp_k_k]
      # puts "-----------------------------"
      if { [catch {drawUncertainty world1 $px $py $ptheta Cp_k_k $scale $points} message] } {
         # puts "Unable to draw $message"
			
         #exit
      } else {
        catch {drawUncertainty world1 $px $py $ptheta Cp_k_k $scale $points}
      }
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
      motionCovariance pCdelta $px $py $ptheta $vl $vr $dt $kl $kr
      #motionJacobian Fdelta $x $y $theta $vl $vr $dt
      motionJacobian pFdelta $px $py $ptheta $vl $vr $dt
      poseJacobian pFp $px $py $ptheta $vl $vr $dt
      # update Cp from Cp and Cdelta and the Jacobians
      # (error propagation law)
      matset Cp_k1_k [matexpr pFp * Cp_k_k * ~pFp + pFdelta * pCdelta * ~pFdelta]
       # get transformation matrix for current pose
      robosim getTransformationMatrix $px $py $ptheta H
      # update pose
      lassign [updatePose $px $py $ptheta $vl $vr $dt] px1 py1 ptheta1
      # get transformation matrix for next pose
      robosim getTransformationMatrix $px1 $py1 $ptheta1 H1


      # # puts [formatMatrix Cp]
      # # puts "\{[mat getRowByRow Cp]\}"
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
          # puts "collision at current position"
          break
      }
      # check collision for predicted position and movement
      set c1 [robosim collision H1]
      set mc [robosim movementCollision H H1]
      if {$c1 || $mc} {
          # puts "predicted collision at $x $y $theta $vl $vr"
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
    ## puts $theta
    set theta_norm [normAngle $theta]
    # puts "Theta $theta_norm"


    mat set p_k1_k 0 0 $px1
    mat set p_k1_k 1 0 $py1
	set ptheta1 [normAngle $ptheta1]
    mat set p_k1_k 2 0 $ptheta1

      #############################################################################
      # 2. Read sensordata from one arbitrary beacon (calculation of alpha=z_tilde_P_k1 and beta=z_tilde_K_k1)
      #############################################################################
      # Get random beacon
    set K [ expr int( rand() * $numBeacons ) ]
    # puts "K $K"
    set atanTmp [expr atan2([lindex $beacon_k($K) 1] - $y, [lindex $beacon_k($K) 0] - $x)]
    set alpha [ expr $atanTmp - $theta_norm - [ gsl randGaussian $sigmaP ] ]
	set alpha [normAngle $alpha]
    # puts "alpha $alpha"
    ## puts "beacon: $beacon_k($K)"
    ## puts "becon: $atanTmp"
    ## puts "compas: $thetaCompas"
    ## puts " alpha_K $alpha_K"
    ## puts $alpha_K

    # Get compas
    set beta [ expr $theta + [ gsl randGaussian $sigmaK ] ]
	set beta [normAngle $beta]
    # puts "beta $beta"
    # # puts "Compas: $beta"
    # # puts "Compas + Noise: $betaNoise"

      #############################################################################
      # 3. Calculate the hypotheses for every beacon i (Calculate z_hat_k1)
      #############################################################################

      # Iterate through every beacon
      set matches 0
	  set current_Gate_value inv
      for {set i 0} {$i < $numBeacons} {incr i} {
          # Get the hypothesis
          set z_hat_k1 [ expr atan2([lindex $beacon_k($i) 1] - $py1, [lindex $beacon_k($i) 0] - $px1) - $beta]
		  set z_hat_k1 [normAngle $z_hat_k1]
          # puts "z_hat_k1 $z_hat_k1"

		  #############################################################################
		  # 4. Calculate innovation
		  #############################################################################
		  # Get Jacobi-Matrix H for the messurement i at current (not groundtruth) position
		  set bx [expr [lindex $beacon_k($i) 0] - $px1]
		  set by [expr [lindex $beacon_k($i) 1] - $py1]
		  set denom [expr 1.0 / ( pow($bx , 2) + pow($by , 2))]

		  set h1 [expr $by / $denom]
		  set h2 [expr $bx / $denom]
		  mat matrix Hi 2 3
		  mat set Hi 0 0 $h1
		  mat set Hi 0 1 $h2
		  mat set Hi 0 2 -1.0
		  mat set Hi 1 0 0.0
		  mat set Hi 1 1 0.0
		  mat set Hi 1 2 1.0



		  # Get covariance R of messurement i
		  mat matrix R 2 2
		  mat set R 0 0 [expr pow($sigmaP,2)]
		  mat set R 0 1 0.0
		  mat set R 1 0 0.0
		  mat set R 1 1 [expr pow($sigmaK,2)]

		  # Get covariance Cv of innovation
		  matset Cv [matexpr Hi * Cp_k1_k * ~Hi + R]

		  # Invert Cv
		  set Cv_rows [mat getRowByRow Cv]
		  set row0 [lindex $Cv_rows 0]
		  set row1 [lindex $Cv_rows 1]
		  set Cv_00 [lindex $row0 0]
		  set Cv_01 [lindex $row0 1]
		  set Cv_10 [lindex $row1 0]
		  set Cv_11 [lindex $row1 1]

		  set det [expr $Cv_00 * $Cv_11 - $Cv_01 * $Cv_10 ]
		  mat matrix Cv_inv 2 2
		  mat set Cv_inv 0 0 [expr $Cv_11 / $det]
		  mat set Cv_inv 0 1 [expr -$Cv_01 / $det]
		  mat set Cv_inv 1 0 [expr -$Cv_10 / $det]
		  mat set Cv_inv 1 1 [expr $Cv_00 / $det]


		  # Get innovation v (call it "nu" to getting confused by velocity)
		  mat matrix nu 2 1
		  mat set nu 0 0 [expr $alpha - $z_hat_k1]
		  mat set nu 1 0 [expr $beta - $ptheta1]
		  # puts [formatMatrix nu]


		  # Calculate the innovation gate and the Kalman gain K
		  set g 0.1
		  matset Gate [matexpr ~nu * Cv_inv * nu]

		  set Gate_value [mat getVectorList Gate]
		  puts "Gate [lindex $Gate_value 0] <= g [expr pow($g,2)]"
		  
		  # if { [lindex $Gate_value 0] <= [expr pow($g,2)] &&  [lindex $Gate_value 0] <= $current_Gate_value} 
		  if { [lindex $Gate_value 0] <= [expr pow($g,2)] } {
		  

			set current_Gate_value $Gate_value
			if {$matches == 0} {
			  incr matches

			  matset K_gain [matexpr Cp_k1_k * ~Hi * Cv_inv]
			  # # puts [formatMatrix Hi]

			 ## puts [formatMatrix Cp_k1_k]
			 ## puts [formatMatrix Hi]
			 ## puts [formatMatrix Cv_inv]
			 ## puts [formatMatrix K_gain]

			} else {
			  incr matches
			}
		  }
      }
      #############################################################################
      # 5. Calculate the new position and variance
      #############################################################################

	  matset Cp_k_k_tmp [matexpr Cp_k1_k - K_gain * Cv * ~K_gain ]
      puts "MATCHES $matches"
      if { $matches == 1 && [checkPosEigVal Cp_k_k]} {
	  # if { $matches == 1} 
         # # puts "Match found"
         ## puts [formatMatrix K_gain]
         # puts "++++++++++++++++++++++++++"
         ## puts [formatMatrix p_k1_k]
         # puts "K_gain: \n[formatMatrix K_gain]"
         ## puts [formatMatrix nu]
         ## puts [formatMatrix Cv]
         ## puts [formatMatrix Cp_k1_k]
		 sleep 1
         matset p_k_k [matexpr p_k1_k + K_gain * nu]
		 # matset p_k_k [matexpr p_k1_k]
         matset Cp_k_k [matexpr Cp_k1_k - K_gain * Cv * ~K_gain ]
         # exit
      } else {
    ## puts "To many/No matches found"
         # puts "++++++++++++++++++++++++++"
         matset Cp_k_k [matexpr Cp_k1_k]
         matset p_k_k [matexpr p_k1_k]
      }
     # sleep 1
    }
    mat leave
}

# uncomment this line for manual control on startup:
if {$doManualControl} {
    manualControlLoop
}