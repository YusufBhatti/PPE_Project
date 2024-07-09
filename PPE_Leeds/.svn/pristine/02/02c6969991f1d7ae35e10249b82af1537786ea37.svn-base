#! /usr/bin/env python
#_________________________________________________________________________________________
#
# Luis kornblueh, MPI-M, 2012
#_________________________________________________________________________________________

from mpmath import *

mp.dps = 50     # calculations with precision of 50 significant decimal digits
pdps = 42       # printing precision of 42 significant decimal digits

print 'REAL(wp), PARAMETER :: euler     = ', nstr(exp(1), pdps), '_wp'
print 'REAL(wp), PARAMETER :: log2e     = ', nstr(log(exp(1), 2), pdps), '_wp'
print 'REAL(wp), PARAMETER :: log10e    = ', nstr(log(exp(1), 10), pdps), '_wp'
print 'REAL(wp), PARAMETER :: ln2       = ', nstr(ln(2), pdps), '_wp'
print 'REAL(wp), PARAMETER :: ln10      = ', nstr(ln(10), pdps), '_wp'
print 'REAL(wp), PARAMETER :: pi        = ', nstr(pi, pdps), '_wp'
print 'REAL(wp), PARAMETER :: pi_2      = ', nstr(fdiv(pi, 2), pdps), '_wp'
print 'REAL(wp), PARAMETER :: pi_4      = ', nstr(fdiv(pi, 4), pdps), '_wp'
#>>SF
print 'REAL(wp), PARAMETER :: pi_6      = ', nstr(fdiv(pi, 6), pdps), '_wp'
#<<SF
print 'REAL(wp), PARAMETER :: rpi       = ', nstr(fdiv(1, pi), pdps), '_wp'
print 'REAL(wp), PARAMETER :: rpi_2     = ', nstr(fdiv(2, pi), pdps), '_wp'
print 'REAL(wp), PARAMETER :: rsqrtpi_2 = ', nstr(fdiv(2, sqrt(pi)), pdps), '_wp'
print 'REAL(wp), PARAMETER :: sqrt2     = ', nstr(sqrt(2), pdps), '_wp'
print 'REAL(wp), PARAMETER :: sqrt1_2   = ', nstr(fdiv(sqrt(2), 2), pdps), '_wp'
print 'REAL(wp), PARAMETER :: sqrt3     = ', nstr(sqrt(3), pdps), '_wp'
print 'REAL(wp), PARAMETER :: sqrt1_3   = ', nstr(fdiv(sqrt(3), 3), pdps), '_wp'
print 'REAL(wp), PARAMETER :: cos45     = ', nstr(sqrt(fdiv(1, 2)), pdps), '_wp'
print 'REAL(wp), PARAMETER :: pi_5      = ', nstr(fdiv(pi, 5), pdps), '_wp'
print 'REAL(wp), PARAMETER :: pi2       = ', nstr(fmul(2, pi), pdps), '_wp'
phi0 = fsub(fdiv(pi, 2),fmul(2, acos(fdiv(1, fmul(2, sin(fdiv(pi, 5)))))))
print 'REAL(wp), PARAMETER :: phi0      = ', nstr(phi0, pdps), '_wp'
print 'REAL(wp), PARAMETER :: rad2deg   = ', nstr(fdiv(180, pi), pdps), '_wp'
print 'REAL(wp), PARAMETER :: deg2rad   = ', nstr(fdiv(pi, 180), pdps), '_wp'
print 'REAL(wp), PARAMETER :: one_third = ', nstr(fdiv(1, 3), pdps), '_wp'
