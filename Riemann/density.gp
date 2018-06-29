reset
set xr [0:1]
set yr [0:6.5]
set ticslevel 0
set size square

set grid
set term gif animate delay 20
set output "density.gif"
unset title
set tics font "Arial,20"
set xl 'time' font "Arial,20"
set yl 'density' font "Arial,20"

n0 = 0
nm = 12
dn = 1

load 'C:\cygwin\home\‹âŽŸ˜Y\cwork\astro\cfd\density_gif.plt'
