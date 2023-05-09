# Dynamic Programming for a simple 1D case.
# Author: Mohamed Naveed Gul Mohamed
# email:mohdnaveed96@gmail.com
# Date: Feb 21st 2021
from __future__ import division

from dp_funcs import *
import math as m
import numpy as np
import time

FILE_WRITE = False# False | True

if FILE_WRITE:
    filename = "/home/naveed/Dropbox/Research/Data/CDC21/stochastic_dp_1d_T60000_X50_processNoise_e0.csv"
    file = open(filename,"a")
    file.write('epsilon' + ',' + 'Average Cost' + ',' + 'Cost variance' + ',' + 'Average Time' + '\n' )

def Jx_func(j, J, X):

    #computing numberical gradients
    if j!=0 and j!=(X.size-1):
        Jx = 0.5*((J[j+1] - J[j])/(X[j+1] - X[j]) + (J[j] - J[j-1])/(X[j] - X[j-1]))
        #print("j=", j," Jx:", Jx," X(j+1):", X[j+1], " X(j):", X[j], " X(j-1):", X[j-1])

    elif j==0:
        Jx = (J[j+1] - J[j])/(X[j+1] - X[j])
        #print("j=", j," Jx:", Jx," X(j+1):", X[j+1], " X(j):", X[j])

    elif j==(X.size-1):
        Jx = (J[j] - J[j-1])/(X[j] - X[j-1])
        #print("j=", j," Jx:", Jx, " X(j):", X[j], " X(j-1):", X[j-1])

    return Jx

def Jxx_func(j, J, X):

    #computing numberical hessians

    if j!=0 and j!=(X.size-1):
        Jxx = 0.5*((Jx_func(j+1,J,X) - Jx_func(j,J,X))/(X[j+1] - X[j]) +
               (Jx_func(j,J,X) - Jx_func(j-1,J,X))/(X[j] - X[j-1]))

    elif j==0:
        Jxx = (Jx_func(j+1,J,X) - Jx_func(j,J,X))/(X[j+1] - X[j])

    elif j==(X.size-1):
        Jxx = (Jx_func(j,J,X) - Jx_func(j-1,J,X))/(X[j] - X[j-1])

    return Jxx

def solve_DP_Bellman(X, x0, xT, n, N, epsilon = 0.0):

    dt = 1.0/N

    J = np.zeros((N+1,n)) #Value function
    J_temp = np.zeros(n)
    u = np.zeros((N,n)) #control

    print('DP epsilon =', epsilon)

    #solving DP
    for t in range(N,-1,-1):

        if t == N:

            for j in range(X.size):
                J[t,j] = terminal_cost(X[j], xT)

        else:

            for j in range(X.size):
                for k in range(X.size):
                    #Calculating the control $u$ to take X(j) to X(k)
                    u_temp = find_u(X[j],X[k],dt)
                    J_temp[k] = current_cost(X[j], u_temp, xT, dt) + Exp_J(X[j], u_temp, J[t+1,:], X, dt, epsilon)

                    #print('k = ', k, 'u =', u_temp, ' inc cost = ', current_cost(X[j],find_u(X[j],X[k],dt),dt), ' J(t+1,k) = ', J[t+1,k])

                x_next_idx = np.argmin(J_temp)


                #print('j = ', j, '   J_temp', J_temp)
                #print('current idx = ', j,'  x_next_idx = ', x_next_idx)

                u[t,j] = find_u(X[j],X[x_next_idx],dt)
                J[t,j] = J_temp[x_next_idx]


    #print "J:\n",J
    #print "u:\n",u

    return J,u

def solve_DP_HJB(X, x0, xT, n, N, epsilon = 0.0):

    dt = 1.0/N

    J = np.zeros((N+1,n)) #Value function
    J_temp = np.zeros(n)
    u = np.zeros((N,n)) #control

    print('DP epsilon =', epsilon)

    #solving DP
    for i in range(N,-1,-1):

        if i == N:

            for j in range(X.size):
                J[i,j] = terminal_cost(X[j], xT)

        else:

            for j in range(X.size):

                #print("i:", i, " j", j)
                Jx = Jx_func(j,J[i+1,:],X)
                u[i,j] = - Jx*g(X[j])/R #optimal control

                Jxx = Jxx_func(j,J[i+1,:],X)

                #print('i = ', i, ' j = ', j, ' Jx = ', Jx, ' Jxx = ', Jxx)

                J[i,j] = current_cost(X[j], u[i,j], xT, dt) + J[i+1,j] + Jx*g(X[j])*u[i,j]*dt +\
                            Jx*f(X[j])*dt + 0.5*Jxx*((epsilon*u_max)**2)*dt

                Courant_num = abs(f(X[j])*dt/(X[1] - X[0]))
                #print('C num = ', Courant_num)

                if Courant_num >=1:
                    print('CFL condition violated!')


    #print "J:\n",J
    #print "u:\n",u

    return J,u

def monteCarloSims(X, x0, xT, N, dt, del_t, epsi_range, iters=50):

    N_steps = int(N*dt / del_t) # number of steps the system propagates

    factor_del_t = int(del_t / dt)

    x_sol = np.zeros(N_steps+1)
    x_sol[0] = x0

    #J, u = solve_DP_Bellman(X, x0, xT, n, N, 0.0)
    J, u = solve_DP_HJB(X, x0, xT, n, N, 0.0)

    for epsilon in epsi_range:

        start_exec = time.time()

        #J, u = solve_DP_HJB(X, x0, xT, n, N, epsilon)
        #J, u = solve_DP_Bellman(X, x0, xT, n, N, epsilon)

        cost_vector = np.zeros(iters)

        for iter in range(iters):

            U_opti = np.zeros(N_steps)

            for i in range(N_steps):

                idx = find_nearest(X,x_sol[i])

                U_opti[i] = u[i*factor_del_t,idx]

                x_sol[i+1] = model(x_sol[i],U_opti[i], del_t, epsilon)


            #print "Solution:", x_sol
            #print "U:", U_opti
            cost = calculate_cost(x_sol,U_opti,xT,del_t)
            cost_vector[iter] = cost

        end_exec = time.time()
        time_taken = (end_exec - start_exec)/iters
        print("epsilon:",epsilon , "   Average Cost:",np.mean(cost_vector), "  Cost var:", np.var(cost_vector))

        if FILE_WRITE:
            file.write(str(epsilon) + ',' + str(np.mean(cost_vector)) + ',' + str(np.var(cost_vector)) + ',' + str(time_taken) + '\n')

def sampleTrial(X, x0, xT, N, dt, del_t, epsilon = 0.0):

    J, u = solve_DP_HJB(X, x0, xT, n, N, 0.0)
    #J, u = solve_DP_Bellman(X, x0, xT, n, N, 0.0)

    print('DP solved.')

    N_steps = int(N*dt / del_t) # number of steps the system pro

    factor_del_t = int(del_t / dt)

    #print('N_steps = ', N_steps, ' Factor = ', factor_del_t)

    x_sol = np.zeros(N_steps+1)
    x_sol[0] = x0
    U_opti = np.zeros(N_steps)

    print('Execution epsilon:', epsilon)

    for i in range(N_steps):

        idx = find_nearest(X,x_sol[i])

        U_opti[i] = u[i*factor_del_t, idx]

        x_sol[i+1] = model(x_sol[i],U_opti[i], del_t, epsilon)

    print("Solution:")
    for t, x_temp, u_temp in zip(range(N_steps),x_sol, U_opti):
        print('t =', t, '   X:', x_temp, '   U:',u_temp)
    print('t =', N_steps, '   X:', x_sol[N_steps])

    cost = calculate_cost(x_sol,U_opti,xT,del_t)
    print("cost:", cost)

    cost_to_go = np.zeros(N_steps+1)
    for i in range(N_steps+1):
        j = find_nearest(X,x_sol[i])
        cost_to_go[i] = J[i*factor_del_t,j]

    plot_func(J,X,N,x_sol,cost_to_go,U_opti,dt)

if __name__=='__main__':

    x0 = 1.0 #initial state.
    xT = 4.8 #final state.

    n = 50 #state space discretization size
    N = 60000 #number of time steps
    dt = 1.0/N #dt - time step for DP solution.
    del_t = 0.02 #del_t - time step for model update.

    X = np.linspace(0,5,n) #space discretization

    #np.random.seed(2)

    sampleTrial(X, x0, xT, N, dt, del_t, epsilon = 0.0)
    #monteCarloSims(X, x0, xT, N, dt, del_t, epsi_range = np.linspace(0.0,1.0,11),iters=100)
