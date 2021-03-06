using LinearAlgebra, Plots, ECOS, Convex, DelimitedFiles, Test, PandaMPC

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

span = 1
resample = 0.02
final_z =  [1; 0; 0; 1; 0; 0; 1; 0; 0; 1; 0; 0; 1; 0; 0; 1; 0; 0 ; 1; 0; 0]
initial_x = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
initial_u = [0; 0; 0; 0; 0; 0; 0]
initial_z = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
pesoQ = 1
pesoR = 0.001
m = 7;  n = 21; l = 21
T = 15+1 # Number of timesteps/ Prediction (and Control) Horizon

tfin, zfin, ufin = go7(span, final_z, initial_x, initial_u, initial_z, resample, T, n, l, m, pesoQ, pesoR)

a = (abs.(zfin[:,end]).-final_z).<0.01

@test a[1]
@test a[2]
@test a[3]
@test a[4]
@test a[5]
@test a[6]
@test a[7]
@test a[8]
@test a[9]
@test a[10]
@test a[11]
@test a[12]
@test a[13]
@test a[14]
@test a[15]
@test a[16]
@test a[17]
@test a[18]
@test a[19]
@test a[20]
@test a[21]

