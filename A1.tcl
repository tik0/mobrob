# ===========================================================================
# 
# A1.tcl --
#
# 0.1
# ===========================================================================

puts "([info script])"

robosim clearWorld
robosim clearRobot

# Setup non-driven static wheels
# wheel radius
set r 0.5
# wheel distance from center
set R 1.0
# wheel width
set w 0.3
# create left and right wheel
set leftWheelID [robosim addWheel $r $w $R $PIHALF 0.0 1 0]
set rightWheelID [robosim addWheel $r $w $R [expr -$PIHALF] $PI 1 0]

# Setup driven castor wheel
# wheel radius
set rc [expr $r]
# wheel distance from center
set L 2.0
# wheel width
set wc [expr $w]
# create front wheel
set frontWheelID [robosim addWheel $rc $wc $L 0.0 $PIHALF 1 1]


# create body
set l [expr $R]
set l2 [expr $l / 2.0]
set ml2 [expr -$l2]
set ml [expr -$l]
# robosim addBodyPart x1   y1 x2 y2
# top-right
robosim addBodyPart $L $l [expr $L+$rc] $l2
# top
robosim addBodyPart $ml2 $l $L $l
# right
robosim addBodyPart [expr $L+$rc] $l2 [expr $L+$rc] $ml2
# bottom-right
robosim addBodyPart [expr $L+$rc] $ml2 $L $ml
# bottom
robosim addBodyPart $l2 $ml $L $ml
# bottom-left
robosim addBodyPart $ml2 $ml $ml $ml2
# left
robosim addBodyPart $ml $ml2 $ml $l2
# top-left
robosim addBodyPart $ml $l2 $ml2 $l
# orientation
robosim addBodyPart 0 0 [expr $L+$rc] 0

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

# path integration,
proc updatePose {x y theta v dt beta L} {
    global l
    set cb [expr cos($beta)]
	set sb [expr sin($beta)]
    set ct [expr cos($theta - $dt / 2.0 / $L * $v * $cb)]
    set st [expr sin($theta - $dt / 2.0 / $L * $v * $cb)]
    set dx_dt [expr $v * $sb * $ct ]
    set dy_dt [expr $v * $sb * $st ]
    set dTheta_dt [expr -1 / $L * $v * $cb ]
    set x [expr $x + $dx_dt * $dt ]
    set y [expr $y + $dy_dt * $dt ]
    set theta [expr $theta + $dTheta_dt * $dt ]
    return "$x $y $theta"
}

proc motionJacobian {FdeltaName x y theta v beta dt} {
    global R L
    set cb [expr cos($beta)]
    set sb [expr sin($beta)]
    set arg [expr $theta - $dt / 2.0 / $L * $v *$cb]
    set ct [expr cos($arg)]
    set st [expr sin($arg)]
    ## F_sv
    # First row of F_sv, left-term and right-term
    set lt1 [expr $dt * $sb * $ct]
    set rt1 [expr $v * $dt * $dt / 2.0 / $L * $cb * $sb * $st]
    # Second row of F_sv, left-term and right-term
    set lt2 [expr $dt * $sb * $st]
    set rt2 [expr - $v * $dt * $dt / 2.0 / $L * $cb * $sb * $ct]
    # Third row of F_sv, left-term and right-term
    set t3 [expr -$dt / $L * $cb]
    ## F_sb
    # First row of F_sv, left-term and right-term
    set lt11 [expr $v * $dt * $cb * $ct]
    set rt11 [expr - $v * $dt * $sb * $dt / 2.0 / $L * $v * $sb * $st]
    # First row of F_sv, left-term and right-term
    set lt22 [expr $v * $dt * $cb * $st]
    set rt22 [expr $v * $dt * $sb * $dt / 2.0 / $L * $v * $sb * $ct]
    # Third row of F_sv, left-term and right-term
    set t33 [expr $v * $dt * $sb / $L]
    
    mat setRowByRow $FdeltaName \
      "{[expr $lt1 + $rt1] [expr $lt11 + $rt11]} \
         {[expr $lt2 + $rt2] [expr $lt22 + $rt22]} \
         {[expr $t3]      [expr $t33]}"
}


# pose Jacobian, see Siegwart/Nourbakhsh p.189
proc poseJacobian {FpName x y theta v beta dt} {
    global R L
    set sb [expr sin($beta)]
    set cb [expr cos($beta)]
    set arg [expr $theta - $dt / 2.0 / $L * $v * $cb]
    set m13 [expr -$dt * $v * $sb * sin($arg)]
    set m23 [expr  $dt * $v * $sb * cos($arg)]
    mat setRowByRow $FpName \
      "{1.0 0.0 [expr $m13]} \
         {0.0 1.0 [expr  $m23]} \
         {0.0 0.0 1.0}"
}

# motion covariance, see Siegwart/Nourbakhsh p.188

proc motionCovariance {CdeltaName v dt dbeta ks kb} {
    set s [expr $dt * $v]
    mat setRowByRow $CdeltaName \
      "{[expr $ks * abs($s)] 0.0} \
         {0.0 [expr $kb * abs($dbeta)]}"
}

#----------------------------------------------------------------------------
# robot is controlled manually
#----------------------------------------------------------------------------

# time step
set dt 0.5
# wheel angular
set beta [expr $PIHALF]
# speed state
set v 0.0
# set vl 0.0
# set vr 0.0
# rotatory/translatory speed increase/decrease
set dbeta 0.1
set betaSet 0
# set dvr 0.05
set dvt 0.05
# used for keyboard control
set stopLoop 0
# factors for motion covariance
set ks 0.1
set kb 0.1
# initial variance in each direction
set c00 0.1
set c11 0.1
set c22 0.01

proc keyHandler {canvasName key} {
      global v beta dbeta dvt stopLoop freeze betaSet
#     global vl vr dvr dvt stopLoop freeze

    switch $key {
      "Left" {
          set beta [expr $beta + $dbeta]
		  set betaSet 1
      }
      "Right" {
          set beta [expr $beta - $dbeta]
		  set betaSet 1
      }
      "Up" {
          set v [expr $v + $dvt]
      }
      "Down" {
          set v [expr $v - $dvt]
      }
      "space" {
          set v 0.0
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

proc manualControlLoop {} {
      global x y theta dt v beta dbeta betaSet L rc wc stopLoop ks kb c00 c11 c22 scale points freeze frontWheelID
    
#     source simple.world
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
    while {!$stopLoop} {
      # draw new robot
      deleteRobot world1
	  robosim setSteering $frontWheelID $beta
      drawRobot world1 $x $y $theta
	  update
      # draw uncertainty ellipse
      deleteEllipse world1
      drawUncertainty world1 $x $y $theta Cp $scale $points
      if {$freeze} {
          drawUncertainty world1 $x $y $theta Cp $scale $points \
            1 black frozen
          set freeze 0
      }
      update
      # update uncertainty
	  if {$betaSet} {
		motionCovariance Cdelta $v $dt $dbeta $ks $kb
		} else {
		motionCovariance Cdelta $v $dt 0.0 $ks $kb
		}
      motionJacobian Fdelta $x $y $theta $v $beta $dt
      poseJacobian Fp $x $y $theta $v $beta $dt
      # update Cp from Cp and Cdelta and the Jacobians
      # (error propagation law)
      matset Cp [matexpr Fp * Cp * ~Fp + Fdelta * Cdelta * ~Fdelta]
      # get transformation matrix for current pose
      robosim getTransformationMatrix $x $y $theta H
      # update pose
      lassign [updatePose $x $y $theta $v $dt $beta $L] x1 y1 theta1
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
          set v 0.0
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
    }
    mat leave
}


manualControlLoop