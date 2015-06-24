#
# FILE: pm3d_setting.gpl
#
# set a default settings for viewing PLUTO binary data 
# with Gnuplot using the pm3d style. 
# It also provides the useful macros @fltform and
# @dblform to conveniently plot data in single or double precision.
#
# Data can be read from a single dataset (in which case you can set nvar=<n> 
# to select the variable) or from mutiple files (in which case you need to 
# supply the desired filename to splot)
#
# Last modified: Feb 15, 2014 by A. Mignone (mignone@ph.unito.it) 
#

set autoscale xfixmin
set autoscale xfixmax
set autoscale yfixmin
set autoscale yfixmax

unset pm3d
set title "  "
set xlabel "x"
set ylabel "y"
#set key top left
unset key
set contour
set cntrparam level 20
unset surface 
set view map

# set margins

set lmargin at screen 0.075
set rmargin at screen 0.85
set bmargin at screen 0.1
set tmargin at screen 0.925

# -- Set colormap macros --

set palette defined

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

# find maximum length and fix the aspect ratio

#Lmax = (Lx > Ly ? Lx:Ly)
#set size Lx/Lmax,Ly/Lmax

