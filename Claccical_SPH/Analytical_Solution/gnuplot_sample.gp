reset;
cd "data"
set size square;
set grid;
#load "../make_animation1.plt"

set term png
load "../make_figure.plt"

!mv -f *.gif ../
!mv -f *.png ../
