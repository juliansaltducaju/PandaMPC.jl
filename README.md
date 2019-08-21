# PandaMPC.jl
The goal of this package is to calculate a trajectory to feed it offline to a Franka Emika Panda Robot, whose ultimate goal is to catch a ball with its end-effector.

The trajectory calculation is calculated using Model Predictive Control, which allows to solve a convex optimization problem by minimizing a cost function under the presence of physical and temporal constraints and in order to achieve the goal position of the end-effector.



