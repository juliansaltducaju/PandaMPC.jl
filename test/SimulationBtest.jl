using Test, PandaMPC

using ECOS
using Convex
using Plots

import PandaMPC

#Initial conditions
x0 = -1;
y0 = 0;
z0 = 0;
v0x = 0.9;
v0y = 0.2;
v0z = 4;
cond0 = [x0; y0; z0; v0x; v0y; v0z];
simtimeball = 1; # Simulation time

# Here we get the trajectory of the ball
t, posdrag, posnodrag = ball(cond0, simtimeball)
xplane = 0.4; #distance of grabbing to the ball from the base
xcapture, ycapture, zcapture, tcapture = intercplane(t, posdrag, xplane)

pesoQ = 1 #state penalization weight
pesoR = 0.001 #input penalization weight
m = 1 #number of joints



Q = [0 0 0; 0 1 0; 0 0 1]
R = pesoR
@test (Q,R) == weights(pesoQ, pesoR, m)
@test (xcapture, ycapture, zcapture, tcapture) == intercplane(t, posdrag, xplane)
