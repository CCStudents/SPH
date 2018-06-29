reset
set xr [0:1]
set yr [0:11000]
set ticslevel 0
set size square

set grid
set term gif animate delay 20
set output "pressure.gif"
unset title
set tics font "Arial,20"
set xl 'time' font "Arial,20"
set yl 'pressure' font "Arial,20"

n0 = 0
nm = 12
dn = 1

load 'C:\cygwin\home\‹âŽŸ˜Y\cwork\astro\cfd\pressure_gif.plt'
