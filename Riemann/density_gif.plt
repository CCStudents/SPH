if(exist ("n") ==0 || n < 0 ) n = n0

title(n) = sprintf("t=%d",n)
set title title(n)

filename = sprintf("exact_%d.dat",n)
plot filename using 1:2 w l notitle

n = n + dn
if ( n < nm ) reread
undefine n
