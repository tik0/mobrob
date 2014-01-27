#!/usr/bin/wish
frame .buFrame
button .buDraw -text Draw -command drawCircle
button .buErase -text Erase -command cleanCanvas
button .buExit -text Exit -command exit

pack .buDraw .buErase .buExit -in .buFrame -side top -anchor n
pack .buFrame -in .

canvas .c -width 300 -height 200
pack .c -in .

proc drawCircle { } {
  .c create oval 70 20 230 180 -fill red -outline {} -tags circle
}

proc cleanCanvas { } {
  .c delete circle
}
