if(exist ("n") ==0 || n < 0 ) n = n0

filename = sprintf("exact_%d.dat",n)
plot filename using 1:3 w l notitle

n = n + dn
if ( n < nm ) reread
undefine n
