# ===========================================================================
# 
# differential.tcl --
# manual control with uncertainty visualization / simple obstacle avoidance
# for a robot with differential steering
#
# Ralf Moeller
# 
#    Copyright (C) 2006
#    Computer Engineering Group
#    Faculty of Technology
#    University of Bielefeld
#    www.ti.uni-bielefeld.de
# 
# 1.0 / 29. Mar 06 (rm)
# - from scratch
# 1.1 / 30. Mar 06 (rm)
# 1.2 / 30. Mar 06 (rm)
# 1.3 / 31. Mar 06 (rm)
# 1.4 / 12. Apr 06 (rm)
# - with matexpr
# 1.5 /  6. Mar 07 (rm)
# - with freeze of uncertainty ellipse (Tab)
# 1.6 / 14. Jun 07 (rm)
# 1.7 / 21. Dec 11 (ak)
# - corrected for errors due to new syntax in robosim
#
# ===========================================================================

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
    global vl vr dvr dvt stopLoop freeze

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

proc manualControlLoop {} {
    global x y theta dt vl vr stopLoop kl kr c00 c11 c22 scale points freeze

    source simple.world
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
          puts "predicted collision at $x1 $y1 $theta1"
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
