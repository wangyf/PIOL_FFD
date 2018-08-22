# Input file for PIOL
debug = 0
dt = 0.002
dx = (25, 25, 25)
itio = 50
itstats = 10
mpin = 1
mpout = 1
nn = (841, 841, 1)
np3 = (1, 1, 1)
nt = 3001
tm0 = 0.0

#User customized field (Far-field displacement)
radius = 1e6 #m
rho = 2670. #kg/m^3
vp = 6000. #m/s
vs = 3640. #m/s

#receiver 
t0s = 160. #sec
t0e = 180.
t1s = 0.
t1e = 0.
t2s = 0.
t2e = 0.
t3s = 0.
t3e = 0.
dtheta = 90.
dphi = 360.
# otheta and ophi are optional (must after dtheta and dphi)
otheta = 30
ophi = -160.


fieldio = [
# define a constant field
#('=', 1, 'const', 1.0, (0.0, 0.0, 0.0), (0.0, 0.0, 0.0), 1, [(1, 41, 1), (1, 41, 1), (1, 2, 1), (0, 0, 1)], '-', 2670.0, ['rho']),
# read or write a file name
('=r', 1, 'const', 1.0, (0.0, 0.0, 0.0), (0.0, 0.0, 0.0), 10, [(1, 841, 1), (1, 841, 1), (1, 1, 1), (1, 3001, 1)], 'mr31', 1.0, ['mr31']),
#('=w', 1, 'const', 1.0, (0.0, 0.0, 0.0), (0.0, 0.0, 0.0), 10, [(1, 161, 1), (1, 1, 1), (1, 2, 1), (1, 801, 1)], 'tinti', 1.0, ['tinti']),
# note:
# 1. static field (it == 0)
#         ii(:,4) = (0,0,1)
#    static field (it == nt)
#         ii(:,4) = (nt,nt,1)
# 2. time-variable (it = from 1,nt)
#         ii(:,4) = (1,nt,1)
]
