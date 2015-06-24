# Set a default settings for viewing PLUTO binary data files
# with Gnuplot using the pm3d style. 
# It also provides the useful macros @fltform and
# @dblform to conveniently plot data in single or double precision.
#
# Data can be read from a single dataset (in which case you can set nvar=<n> 
# to select the variable) or from mutiple files (in which case you need to 
# supply the desired filename to splot)
#
# Last Modified: April 22, 2014 by A. Mignone (mignone@ph.unito.it) 
#

print "> Load default pm3d setting"

set autoscale xfixmin
set autoscale xfixmax
set autoscale yfixmin
set autoscale yfixmax

unset contour
set title "  "
set xlabel "x"
set ylabel "y"
#set key top left
unset key
set pm3d map

# -- Set margins --

set lmargin at screen 0.15
set rmargin at screen 0.8
set bmargin at screen 0.1
set tmargin at screen 0.93

# -- Set colormap macros --

set palette defined
set macro
HOT = "rgbformulae 22,13,-31"
RED = "rgbformulae 21,22,23"
RYG = 'model RGB defined ( 0 "red", 0.5 "yellow", 1 "green" )'
RAINBOW = "rgbformulae 33,13,10 " # rainbow (blue-green-yellow-red)
JET = "defined (0  0.0 0.0 0.5, \
                1  0.0 0.0 1.0, \
                2  0.0 0.5 1.0, \
                3  0.0 1.0 1.0, \
                4  0.5 1.0 0.5, \
                5  1.0 1.0 0.0, \
                6  1.0 0.5 0.0, \
                7  1.0 0.0 0.0, \
                8  0.5 0.0 0.0 )"

GREEN_CT = 'model RGB defined (0 "black", 1 "slateblue1", 2 "white")'

# -- Do some printing --

print "  - available colormap macros: HOT, RED, RAINBOW, JET"
