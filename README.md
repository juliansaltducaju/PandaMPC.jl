# PandaMPC.jl

[![Build Status](https://travis-ci.org/juliansaltducaju/PandaMPC.jl.svg?branch=master)](https://travis-ci.org/juliansaltducaju/PandaMPC.jl)


The goal of this package is to calculate a trajectory to feed it offline to a Franka Emika Panda Robot, whose ultimate goal is to catch a ball with its end-effector.

The trajectory calculation is calculated using Model Predictive Control, which allows to solve a convex optimization problem by minimizing a cost function under the presence of physical and temporal constraints and in order to achieve the goal position of the end-effector.

In the following example, the 7 joints of the Panda Robot are going to be driven from an angular postion on 0 radians and no velocity or acceleration to an angular position of 1 radians and also no velocity of acceleration. This will be done within 1 second:

```
m = 7;  n = 21; l = 21
 
initial_x = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]
initial_u = [0; 0; 0; 0; 0; 0; 0]
initial_z = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]

final_z =  [1; 0; 0; 1; 0; 0; 1; 0; 0; 1; 0; 0; 1; 0; 0; 1; 0; 0 ; 1; 0; 0]

span = 1
```
As the trajectory design parameters, we will choose a prediction horizon of 15 steps in the MPC and a penalization on the states 1000 times bigger than the penalization on the inputs:

```
T = 15+1
pesoQ = 1
pesoR = 0.001 
```
The last design parameter will be the `resample`, which allows us to increase the accuracy of the trajectory by recalculating it:
```
resample = 0.02 #The trajectory will be recalculated every 0.02 seconds
```
Finally, to generate the trajectory:
```
t, z, u = go7(span, final_z, initial_x, initial_u, initial_z, resample, T, n, l, m, pesoQ, pesoR)

```
![Results](/example/fig.png)
