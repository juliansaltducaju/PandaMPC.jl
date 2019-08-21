using Test, JuliansMPCforJulia

using ECOS
using Convex
using Plots

import JuliansMPCforJulia: weights, model, go, iteraciones2, plotjoint

pesoQ = 1 #state penalization weight
pesoR = 0.001 #input penalization weight
m = 1 #number of joints



Q = [0 0 0; 0 1 0; 0 0 1]
R = pesoR
@test (Q,R) == weights(pesoQ, pesoR, m)
