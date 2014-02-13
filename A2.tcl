# ===========================================================================
# 
# A2.tcl --
#
# 0.1
# ===========================================================================
puts "([info script])"

robosim clearWorld
robosim clearRobot

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
    set xR 0.0
    set yR 0.0
    set phiR 0.0
    set sweep $PI
    set nRays 21
    set range 10.0
    set sigmaSlope 0.0
    set sigmaOffset 0.2
    robosim createLaser $xR $yR $phiR $sweep $nRays $range $sigmaSlope $sigmaOffset

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
set x 5.0
set y 5.0
set theta [expr -$PI / 2]

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
# robot is doing a simple obstacle avoidance
#----------------------------------------------------------------------------

# time step
set dt 0.1

set vpos 0.5
set vneg -0.5

# Parameter for wallfollowing
set wallDistance 2.5
set thresholdDist 1.5


proc wallFollowing {} {
    global x y theta dt nRays range thresholdDist vpos vneg stopLoop wallDistance

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
      drawLaserScan world1 $scan
      update
      #####################################
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
      
      # find distance on the left side
      set beam [lindex $scan 0]
      set dist [lindex $beam 4]
      set leftDist $leftMinDist
      if {$dist < $leftMinDist} {
          set leftDist $dist
          puts "Left distance: $leftDist"
      }
      puts "Left distance: $dist"
      
      
      # find closest distance on the right side
      set rightMinDist $range
      for {set n [expr $nRays / 2]} {$n < $nRays} {incr n} {
          set beam [lindex $scan $n]
          set dist [lindex $beam 4]
          if {$dist < $rightMinDist} {
            set rightMinDist $dist
          }
      }
      # find distance on the right side
      set beam [lindex $scan [expr $nRays - 1]]
      set dist [lindex $beam 4]
      set rightDist $rightMinDist
      if {$dist < $rightMinDist} {
        set rightDist $dist
        puts "Right distance: $rightDist"
      }
      puts "Right distance: $dist"
      
      # analyse both sides
      set leftWallClose [expr $leftDist < $wallDistance]
      set rightWallClose [expr $rightDist < $wallDistance]
      set leftWallFar [expr $leftDist > $wallDistance]
      set rightWallFar [expr $rightDist > $wallDistance]
      
      set leftClose [expr $leftMinDist < $thresholdDist]
      set rightClose [expr $rightMinDist < $thresholdDist]
      if {$leftClose && $rightClose} {
          set vl $vpos
          set vr $vneg
      } elseif {$leftClose} {
          set vl $vpos
          set vr $vneg
      } elseif {$rightClose} {
          set vl $vpos
          set vr $vneg
      } else {
          set vl [expr 2 * $vpos]
          set vr [expr 2 * $vpos]
          
          if {$leftWallClose} {
            set vl [expr 1.5 * $vpos]
            set vr $vpos
          } elseif {$leftWallFar} {
            set vl $vpos
            set vr [expr 2.5 * $vpos]
          } 
          
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
    # source "obstacleAvoidance.world"
    # source "simple.world"
    source "TestWorld.world"
    drawWorld world1
    wallFollowing

