module MPC

import LinearAlgebra
import Plots
import ECOS
import Convex
import DelimitedFiles

###############################################################################
# This block contains the ball's trajectory simulation functions
"""
Given an array with the initial position and velocity of the ball [x;y;z;vx;vy;vz],\n
this function calculates the trajectory of the ball (with and without drag force)\n
during the simulation time defined by the user
"""
function ball(cond0, simtime)
    x0 = cond0[1];
    y0 = cond0[2];
    z0 = cond0[3];
    v0x = cond0[4];
    v0y = cond0[5];
    v0z = cond0[6];
    #Required ball data
    m = 0.05; # mass of the ball [kg]
    r = 0.05; # radius of the ball [meters]
    g = 9.81; # gravity
    t = range(0,stop=simtime,length=100)

    C = 0.96211*r^2/m;
    x = zeros(length(t),1); # x initialization
    y = zeros(length(t),1); # y initialization
    z = zeros(length(t),1); # y initialization
    xnodrag = zeros(length(t),1); # x initialization
    ynodrag = zeros(length(t),1); # y initialization
    znodrag = zeros(length(t),1); # y initialization

    #ODE coefficients
    # y direction (vertical)
    c1 = atan((-v0z*sqrt(C*g))/g)/sqrt(C*g);
    c2 = z0 - log(cos(c1*sqrt(C*g)))/C;
    # x direction
    c3 = 1/v0x;
    c4 = x0 - log(c3)/C;
    # z direction
    c5 = 1/v0y;
    c6 = y0 - log(c5)/C;

    for i = 1:length(t)
        x[i] = log(c3 + C*t[i])/C + c4;
        y[i] = log(c5 + C*t[i])/C + c6;
        z[i] = log(cos((c1+t[i])*sqrt(C*g)))/C + c2;

        xnodrag[i] = x0 + v0x*t[i];
        ynodrag[i] = y0 + v0y*t[i];
        znodrag[i] = z0 + v0z*t[i] - 0.5*g*t[i]^2;
    end

    posdrag = [x y z]
    posnodrag = [xnodrag ynodrag znodrag]
    return t, posdrag, posnodrag
end
"""
Given an array with the trajectory of the ball, this function calculates the\n
time remaining until the ball crosses a plane perpendicular to the x-direction\n
of the base frame of the robot. This plane will be at a distance xplane to the\n
origin of the frame. Cartesian position of this point is also provided
"""
function intercplane(t, posdrag, xplane = 0.4)
    xcapture = 0;
    ycapture = 0;
    zcapture = 0;
    tcapture = 0;
    for i in 1:length(t)
        if abs(posdrag[i,1]+xplane) < 0.01
            xcapture = -posdrag[i,1];
            ycapture = posdrag[i,2];
            zcapture = posdrag[i,3];
            tcapture = t[i];
        end
    end
    return xcapture, ycapture, zcapture, tcapture
end






################################################################################
# This block contains the inverse kinematics functions for the Franka Panda Robot
"""
Calculates the inverse kinematics of the Panda Robot. Providing the desired\n
cartesian position of the end-effector and the starting configuration, the function\n
will calculate the angles of each of the 7 joints that will achieve this.
"""
function getjointangles(posicCartes, theta0)
    theta = theta0;
    xd = posicCartes;
    M = [1 0 0 0.088; 0 -1 0 0; 0 0 -1 0.8226; 0 0 0 1];
    S = [0 0 0 0 0 0 0; 0 1 0 -1 0 -1 0; 1 0 1 0 1 0 -1; 0 -0.333 0 0.649 0 1.033 0; 0 0 0 0 0 0 0.088; 0 0 0 -0.0825 0 0 0];

    k = 1;
    while k<800
        res =  zeros(4,4,length(theta))
        for i in 1 : length(theta)
            Sw = S[1:3, i];
            Sv = S[4:6, i];
            Smat = [0 -Sw[3] Sw[2] Sv[1]; Sw[3] 0 -Sw[1] Sv[2]; -Sw[2] Sw[1] 0 Sv[3]; 0 0 0 0];
            res[:,:,i] = exp(Smat*theta[i]);
        end


        T = res[:,:,1]*res[:,:,2]*res[:,:,3]*res[:,:,4]*res[:,:,5]*res[:,:,6]*res[:,:,7]*M;

        Js = zeros(6,7)
        Adtotal = zeros(4,4,6)
        # Joint 1
        Js[:,1] = S[:,1];

        # All the other joints
        Adtotal[:,:,1] = res[:,:,1];
        Adtotal[:,:,2] = res[:,:,1]*res[:,:,2];
        Adtotal[:,:,3] = res[:,:,1]*res[:,:,2]*res[:,:,3];
        Adtotal[:,:,4] = res[:,:,1]*res[:,:,2]*res[:,:,3]*res[:,:,4];
        Adtotal[:,:,5] = res[:,:,1]*res[:,:,2]*res[:,:,3]*res[:,:,4]*res[:,:,5];
        Adtotal[:,:,6] = res[:,:,1]*res[:,:,2]*res[:,:,3]*res[:,:,4]*res[:,:,5]*res[:,:,6];

        for j = 1:size(Adtotal,3)
            Ad = Adtotal[:,:,j];
            R = Ad[1:3,1:3];
            p  = Ad[1:3,4];
            p2 = [0 -p[3] p[2]; p[3] 0 -p[1]; -p[2] p[1] 0];
            Js[:,j+1] = [R zeros(3,3); p2*R R]*S[:,j+1];
        end

        Jss = [Js[4:6,:]; Js[1:3,:]];
        Jbb = adi(T)*Jss;
        Jsb = [T[1:3,1:3] zeros(3,3); zeros(3,3) T[1:3,1:3]] * Jbb; #La buena


        e = [xd[1:3]-T[1:3,4]; 0;0;0];
        theta = theta + pinv(Jsb)*e;
        if abs(e[1])+abs(e[2])+abs(e[3]) < 0.005
            k = 5000;
        end
        k = k+1;
    end

    for i in 1:length(theta)
        theta[i] = theta[i] % (2*pi); #para evitar dar varias vueltas al eje (no afecta a T)
    end
    return theta
end

"""Calculates the adjoint of the inverse of a tranformation matrix T, which is also the inverse of the adjoint of T"""
function adi(T)
    # Computes the adjoint of transformation matrix inv(T)
    R = T[1:3,1:3]'
    t = skew(-R*T[1:3,4])
    Z = zeros(3,3)
    A = [R t*R; Z R]
    return A
end

skew(s)::AbstractMatrix{eltype(s)} = [0 -s[3] s[2];s[3] 0 -s[1]; -s[2] s[1] 0]

"""
Calculates the forward kinematics of the Panda Robot. An array with the 7 joint\n
angles must be provided so that the function calculates the cartesian position\n
of the end-effector
"""
function forwardkin(theta)
    M = [1 0 0 0.088; 0 -1 0 0; 0 0 -1 0.8226; 0 0 0 1];
    S = [0 0 0 0 0 0 0; 0 1 0 -1 0 -1 0; 1 0 1 0 1 0 -1; 0 -0.333 0 0.649 0 1.033 0; 0 0 0 0 0 0 0.088; 0 0 0 -0.0825 0 0 0];

    k = 1;

    res =  zeros(4,4,length(theta))
    for i in 1 : length(theta)
        Sw = S[1:3, i];
        Sv = S[4:6, i];
        Smat = [0 -Sw[3] Sw[2] Sv[1]; Sw[3] 0 -Sw[1] Sv[2]; -Sw[2] Sw[1] 0 Sv[3]; 0 0 0 0];
        res[:,:,i] = exp(Smat*theta[i]);
    end
    T = res[:,:,1]*res[:,:,2]*res[:,:,3]*res[:,:,4]*res[:,:,5]*res[:,:,6]*res[:,:,7]*M;

    return T
end







################################################################################
# This block contains the MPC functions
"""
Generates the matrixes Q and R for solving the problem: quadform(z,Q) + quadform(u,R)\n
    m is the number of joints considered in the robot (1, 2 or 3)
"""
function weights(weightQ, weightR, m)
    if m == 1
        Q = zeros(3,3)
        for i = [2,3]
            Q[i,i] = weightQ
        end
        R = weightR
    elseif m == 2
        Q = zeros(6,6)
        for i = [2,3,5,6]
            Q[i,i] = weightQ
        end
        R = zeros(2,2)
        for j = 1:2
            R[j,j] = weightR
        end
    elseif m == 3
        Q = zeros(9,9)
        for i = [2,3,5,6,8,9]
            Q[i,i] = weightQ
        end
        R = zeros(3,3)
        for j = 1:3
            R[j,j] = weightR
        end
        return Q, R
    elseif m == 7
        Q = zeros(21,21)
        for i = [2,3,5,6,8,9,11,12,14,15,17,18,20,21]
            Q[i,i] = weightQ
        end
        R = zeros(7,7)
        for j = 1:7
            R[j,j] = weightR
        end
        return Q, R
    else
        error("Number of joints should be 1, 2, 3 or 7")
    end
    return Q, R
end

"""
 Generates the model matrixes A, B and Dz\n
    h is the sampling period\n
    n is the number of degrees of freedom (3 per joint; 1, 2 or 3 joints)
"""
function model(h, n)
    Fitilde = [1 h h^2/2; 0 1 h; 0 0 1]
    Gamma1tilde = [h^4/24; h^3/6; h^2/2];
    Gammatilde = [h^3/6; h^2; h];
    if n == 3
        Fi = Fitilde
        Gamma1 = Gamma1tilde
        Gamma = Gammatilde
    elseif n == 6
        Fi = [Fitilde zeros(3,3); zeros(3,3) Fitilde]
        Gamma1 = [Gamma1tilde zeros(3,1); zeros(3,1) Gamma1tilde]
        Gamma = [Gammatilde zeros(3,1); zeros(3,1) Gammatilde]
    elseif n == 9
        Fi = [Fitilde zeros(3,3) zeros(3,3); zeros(3,3) Fitilde zeros(3,3); zeros(3,3) zeros(3,3) Fitilde]
        Gamma1 = [Gamma1tilde zeros(3,1) zeros(3,1); zeros(3,1) Gamma1tilde zeros(3,1); zeros(3,1) zeros(3,1) Gamma1tilde]
        Gamma = [Gammatilde zeros(3,1) zeros(3,1); zeros(3,1) Gammatilde zeros(3,1); zeros(3,1) zeros(3,1) Gammatilde]
    elseif n == 21
        Fi = [Fitilde zeros(3,3) zeros(3,3) zeros(3,3) zeros(3,3) zeros(3,3) zeros(3,3);
        zeros(3,3) Fitilde zeros(3,3) zeros(3,3) zeros(3,3) zeros(3,3) zeros(3,3);
        zeros(3,3) zeros(3,3) Fitilde zeros(3,3) zeros(3,3) zeros(3,3) zeros(3,3);
        zeros(3,3) zeros(3,3) zeros(3,3) Fitilde zeros(3,3) zeros(3,3) zeros(3,3);
        zeros(3,3) zeros(3,3) zeros(3,3) zeros(3,3) Fitilde zeros(3,3) zeros(3,3);
        zeros(3,3) zeros(3,3) zeros(3,3) zeros(3,3) zeros(3,3) Fitilde zeros(3,3);
        zeros(3,3) zeros(3,3) zeros(3,3) zeros(3,3) zeros(3,3) zeros(3,3) Fitilde]

        Gamma1 = [Gamma1tilde zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1);
        zeros(3,1) Gamma1tilde zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1);
        zeros(3,1) zeros(3,1) Gamma1tilde zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1);
        zeros(3,1) zeros(3,1) zeros(3,1) Gamma1tilde zeros(3,1) zeros(3,1) zeros(3,1);
        zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) Gamma1tilde zeros(3,1) zeros(3,1);
        zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) Gamma1tilde zeros(3,1);
        zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) Gamma1tilde]

        Gamma = [Gammatilde zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1);
        zeros(3,1) Gammatilde zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1);
        zeros(3,1) zeros(3,1) Gammatilde zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1);
        zeros(3,1) zeros(3,1) zeros(3,1) Gammatilde zeros(3,1) zeros(3,1) zeros(3,1);
        zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) Gammatilde zeros(3,1) zeros(3,1);
        zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) Gammatilde zeros(3,1);
        zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1) Gammatilde]
    else
        error("Number of joints should be 1, 2, 3 or 7")
    end

    eye9 = zeros(n,n)
    for i = 1:n
        eye9[i,i] = 1
    end

    A = Fi
    B = Gamma + (Fi - eye9)*Gamma1/h
    Dz = Gamma1/h

    return A,B,Dz
end

"""
Robot Timescale\n
    You introduce a list of velocities (and times) with a time spacing greater than 1 milisecond\n
     and this function returns a list of velocities with a exact spacing of 1 milisecond\n
     (the gaps between the original time samples are filled by means of a linear interpolation) \n
"""
function robot_timescale(tfin, xfin2)
    for k = 1:length(tfin)
        tfin[k] = round(tfin[k]; digits=3)
        xfin2[k] = round(xfin2[k]; digits=5)
    end

    trobot = collect(tfin[1]:0.001:tfin[end])
    vrobot = [0;1]

    for i = 1:length(tfin)-1
        samples = tfin[i]:0.001:tfin[i+1]
        for j = 1:length(samples)-1
            vrobot = [vrobot; ((xfin2[i+1]-xfin2[i])/(tfin[i+1]-tfin[i]))*(samples[j] - tfin[i]) + xfin2[i]]
        end
    end
    vrobot = [vrobot[3:end]; xfin2[end]]

    return trobot, vrobot
end

"""
MPC generation\n
    span is the calculation time\n
    resample lets you set how often does the MPC recalculate\n
    final_z is the final state\n
    initial_x, initial_u, initial_z are the initial conditions\n\n

    Example - 2 joints, 6 degrees of freedom and time span of 1 second and recalculated every 20ms:\n
    span = 1\n
    resample = 0.02\n
    final_z =  [1; 0; 0; 1; 0; 0]\n
    initial_x = [0; 0; 0; 0; 0; 0]\n
    initial_u = [0; 0]\n
    initial_z = [0; 0; 0; 0; 0; 0]
"""
function go(span, final_z, initial_x, initial_u, initial_z, resample)
    recalc_times = collect(0:resample:span-resample)
    tfin = [0;1] #only for creating these variables
    xfin = [0 1; 0 1; 0 1; 0 1; 0 1; 0 1]
    ufin = [0 1;0 1]
    for k = 1:length(recalc_times)

        tsim = span - recalc_times[k] #Simulation time
        h = tsim/(T-1)
        (A, B, Dz) = model(h, n)

        # Start and end point
        if k == 1
            t = 0
            for i = 1:(T-1)
                t = [t; t[end]+h]
            end
        else
            x = nothing
            u = nothing
            z = nothing
        end

        # Create a optimization variables of size n x 1.
        x = Variable(n, T)
        u = Variable(m, T)
        z = Variable(l, T)

        # Create a problem instance
        (Q,R) = weights(weightQ, weightR, m);
        problem = minimize(quadform(z,Q) + quadform(u,R))

        # Equality constraints
        problem.constraints += x[:,1] == initial_x
        problem.constraints += u[:,1] == initial_u
        problem.constraints += z[:,1] == initial_z
        problem.constraints += z[:,T] == final_z

        for i in 1 : T -1
          problem.constraints += u[:,i+1] <= [1500; 1500]
          problem.constraints += u[:,i+1] >= [-1500; -1500]
          problem.constraints += z[:,i+1] <= [2; 3.14159; 45; 2; 3.14159; 45]
          problem.constraints += z[:,i+1] >= [-2; -3.14159; -45; -2; -3.14159; -45]
          problem.constraints += x[:,i+1] == A*x[:, i] + B* u[:, i]
          problem.constraints += z[:,i+1] == x[:,i+1] + Dz* u[:,i+1]
        end

        # Solve the problem by calling solve!
        solve!(problem, ECOSSolver())

        # Check the status of the problem
       problem.status # :Optimal, :Infeasible, :Unbounded etc.

        # Get the optimum value
        problem.optval

        # global itera, tbien
        itera = 0
        t3 = collect(recalc_times[k]:h:span)
        for i = 1:length(t3)
            if t3[i] <= resample*k
                tfin = [tfin; t3[i]]
                xfin = [xfin x.value[:,i]]
                ufin = [ufin u.value[:,i]]
                itera = i
            end
        end
        tbien = t3[1:itera]
        if tbien[end] == resample*k
            initial_x = x.value[:,itera]
            initial_u = u.value[:,itera]
            initial_z = z.value[:,itera]
        else #Lineal interpolation
            initial_x = ((x.value[:,itera+1]-x.value[:,itera])/(t3[itera+1]-t3[itera]))*(resample*k-t3[itera]) + x.value[:,itera]
            initial_u = ((u.value[:,itera+1]-u.value[:,itera])/(t3[itera+1]-t3[itera]))*(resample*k-t3[itera]) + u.value[:,itera]
            initial_z = ((z.value[:,itera+1]-z.value[:,itera])/(t3[itera+1]-t3[itera]))*(resample*k-t3[itera]) + z.value[:,itera]
        end
    end

    tfin = tfin[3:end]
    xfin = xfin[:,3:end]
    ufin = ufin[:,3:end]
    return tfin, xfin, ufin
end



"""
Simulation for 2 joints:\n\n

The idea behind this is that the robot first waits for an instruction,\n
then it executes that, waits in that position, goes back to the origin\n
and then waits for the next order.\n\n

pos1 is the array of positions of the first joint (where its first value is the starting point)\n

pos2 is the array of positions of the second joint (where its first value is the starting point)\n

span0 is the initial waiting time

span1 is the time span for getting to the goal state

wait_time is the waiting time at the goal state

span2 is the time span for getting back to the goal state from to the goal state

total_time is the total time of the iteration and in this context determines\n
the waiting time until the next iteration

"""
function iteraciones2(pos1, pos2, span0, span1, wait_time, span2, total_time)
    #Initial Pause - Awaiting instructions
    #span0 = 1
    final_z =  [pos1[1]; 0; 0; pos2[1]; 0; 0]
    initial_x = [pos1[1]; 0; 0; pos2[1]; 0; 0]
    initial_u = [0; 0]
    initial_z = [pos1[1]; 0; 0; pos2[1]; 0; 0]
    (tfin, xfin, ufin) = go(span0, final_z,initial_x, initial_u, initial_z, resample)
    timecalc = tfin
    x = xfin
    u = ufin
    for i = 2:length(pos1)
        #Instruction
        final_z =  [pos1[i]; 0; 0; pos2[i]; 0; 0]
        initial_x = [pos1[1]; 0; 0; pos2[1]; 0; 0]
        initial_u = [0; 0]
        initial_z = [pos1[1]; 0; 0; pos2[1]; 0; 0]
        (tfin, xfin, ufin) = go(span1, final_z,initial_x, initial_u, initial_z, resample)
        for i = 1:length(tfin)
            tfin[i] = tfin[i] + timecalc[end] + 0.0001
        end
        timecalc = [timecalc; tfin]
        x = [x xfin]
        u = [u ufin]

        #Waiting time in goal
        final_z =  [pos1[i]; 0; 0; pos2[i]; 0; 0]
        initial_x = [pos1[i]; 0; 0; pos2[i]; 0; 0]
        initial_u = [0; 0]
        initial_z = [pos1[i]; 0; 0; pos2[i]; 0; 0]
        (tfin, xfin, ufin) = go(wait_time, final_z,initial_x, initial_u, initial_z, resample)
        for i = 1:length(tfin)
            tfin[i] = tfin[i] + timecalc[end] + 0.0001
        end
        timecalc = [timecalc; tfin]
        x = [x xfin]
        u = [u ufin]

        #Back to origin
        final_z =  [pos1[1]; 0; 0; pos2[1]; 0; 0]
        initial_x = [pos1[i]; 0; 0; pos2[i]; 0; 0]
        initial_u = [0; 0]
        initial_z = [pos1[i]; 0; 0; pos2[i]; 0; 0]
        (tfin, xfin, ufin) = go(span2, final_z,initial_x, initial_u, initial_z, resample)
        for i = 1:length(tfin)
            tfin[i] = tfin[i] + timecalc[end] + 0.0001
        end
        timecalc = [timecalc; tfin]
        x = [x xfin]
        u = [u ufin]

        #Waiting next instruction
        span3 = total_time - span1 - span2
        final_z =  [pos1[1]; 0; 0; pos2[1]; 0; 0]
        initial_x = [pos1[1]; 0; 0; pos2[1]; 0; 0]
        initial_u = [0; 0]
        initial_z = [pos1[1]; 0; 0; pos2[1]; 0; 0]
        (tfin, xfin, ufin) = go(span3, final_z,initial_x, initial_u, initial_z, resample)
        for i = 1:length(tfin)
            tfin[i] = tfin[i] + timecalc[end] + 0.0001
        end
        timecalc = [timecalc; tfin]
        x = [x xfin]
        u = [u ufin]
    end
    return timecalc, x, u
end


"""
Plotting for position, velocity, acceleration and input (jerk) of one of the two controlled joints\n
    joint should be 1 or 2

"""
function plotjoint(timecalc, x, joint)
    if joint == 1
        plot(timecalc,[x[1,:],x[2,:],x[3,:],u[1,:]], layout = 4, legend = false)
    elseif joint == 2
        plot(timecalc,[x[4,:],x[5,:],x[6,:],u[2,:]], layout = 4, legend = false)
    elseif joint == 3
        plot(timecalc,[x[7,:],x[8,:],x[9,:],u[3,:]], layout = 4, legend = false)
    elseif joint == 4
        plot(timecalc,[x[10,:],x[11,:],x[12,:],u[4,:]], layout = 4, legend = false)
    elseif joint == 5
        plot(timecalc,[x[13,:],x[14,:],x[15,:],u[5,:]], layout = 4, legend = false)
    elseif joint == 6
        plot(timecalc,[x[16,:],x[17,:],x[18,:],u[6,:]], layout = 4, legend = false)
    elseif joint == 7
        plot(timecalc,[x[19,:],x[20,:],x[21,:],u[7,:]], layout = 4, legend = false)
    else
        error("Choose a valid joint")
    end
end

"""
MPC generation\n
    span is the calculation time\n
    resample lets you set how often does the MPC recalculate\n
    final_z is the final state\n
    initial_x, initial_u, initial_z are the initial conditions\n\n

    Example - 2 joints, 6 degrees of freedom and time span of 1 second and recalculated every 20ms:\n
    span = 1\n
    resample = 0.02\n
    final_z =  [1; 0; 0; 1; 0; 0]\n
    initial_x = [0; 0; 0; 0; 0; 0]\n
    initial_u = [0; 0]\n
    initial_z = [0; 0; 0; 0; 0; 0]
"""
function go7(span, final_z, initial_x, initial_u, initial_z, resample)
    recalc_times = collect(0:resample:span-resample)
    tfin = [0;1] #only for creating these variables
    xfin = [0 1; 0 1; 0 1; 0 1; 0 1; 0 1; 0 1; 0 1; 0 1; 0 1; 0 1; 0 1; 0 1; 0 1; 0 1 ; 0 1; 0 1; 0 1; 0 1; 0 1; 0 1]
    ufin = [0 1; 0 1; 0 1; 0 1; 0 1; 0 1; 0 1]
    for k = 1:length(recalc_times)

        tsim = span - recalc_times[k] #Simulation time
        h = tsim/(T-1)
        (A, B, Dz) = model(h, n)

        # Start and end point
        if k == 1
            t = 0
            for i = 1:(T-1)
                t = [t; t[end]+h]
            end
        else
            x = nothing
            u = nothing
            z = nothing
        end

        # Create a optimization variables of size n x 1.
        x = Variable(n, T)
        u = Variable(m, T)
        z = Variable(l, T)

        # Create a problem instance
        (Q,R) = weights(weightQ, weightR, m);
        problem = minimize(quadform(z,Q) + quadform(u,R))

        # Equality constraints
        problem.constraints += x[:,1] == initial_x
        problem.constraints += u[:,1] == initial_u
        problem.constraints += z[:,1] == initial_z
        problem.constraints += z[:,T] == final_z

        for i in 1 : T -1
          problem.constraints += u[:,i+1] <= [5000; 3000; 4000; 4000; 5000; 7000; 7000]
          problem.constraints += u[:,i+1] >= [-5000; -3000; -4000; -4000; -5000; -7000; -7000]
          problem.constraints += z[:,i+1] <= [2.5; 2; 13; 1.5; 2; 6; 2.5; 2; 8; -0.069; 2; 10; 2.5; 2.5; 13; 3.5; 2.5; 18; 2.5; 2.5; 18]
          problem.constraints += z[:,i+1] >= [-2.5; -2; -13; -1.5; -2; -6; -2.5; -2; -8; -3; -2; -10; -2.5; -2.5; -13; 0; -2.5; -18; -2.5; -2.5; -18]
          problem.constraints += x[:,i+1] == A*x[:, i] + B* u[:, i]
          problem.constraints += z[:,i+1] == x[:,i+1] + Dz* u[:,i+1]
        end

        # Solve the problem by calling solve!
        solve!(problem, ECOSSolver())

        # Check the status of the problem
       problem.status # :Optimal, :Infeasible, :Unbounded etc.

        # Get the optimum value
        problem.optval

        # global itera, tbien
        itera = 0
        t3 = collect(recalc_times[k]:h:span)
        for i = 1:length(t3)
            if t3[i] <= resample*k
                tfin = [tfin; t3[i]]
                xfin = [xfin x.value[:,i]]
                ufin = [ufin u.value[:,i]]
                itera = i
            end
        end
        tbien = t3[1:itera]
        if tbien[end] == resample*k
            initial_x = x.value[:,itera]
            initial_u = u.value[:,itera]
            initial_z = z.value[:,itera]
        else #Lineal interpolation
            initial_x = ((x.value[:,itera+1]-x.value[:,itera])/(t3[itera+1]-t3[itera]))*(resample*k-t3[itera]) + x.value[:,itera]
            initial_u = ((u.value[:,itera+1]-u.value[:,itera])/(t3[itera+1]-t3[itera]))*(resample*k-t3[itera]) + u.value[:,itera]
            initial_z = ((z.value[:,itera+1]-z.value[:,itera])/(t3[itera+1]-t3[itera]))*(resample*k-t3[itera]) + z.value[:,itera]
        end
    end

    tfin = tfin[3:end]
    xfin = xfin[:,3:end]
    ufin = ufin[:,3:end]
    return tfin, xfin, ufin
end

end
