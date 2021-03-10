# Dynamic Programming for a simple 1D case.
# Author: Mohamed Naveed Gul Mohamed
# email:mohdnaveed96@gmail.com
# Date: Feb 21st 2021
from __future__ import division

from dp_funcs import *
import math as m
import numpy as np


FILE_WRITE = False# False | True

if FILE_WRITE:
    filename = "/home/naveed/Documents/Dynamic_programming/stochastic_dp_1d.csv"
    file = open(filename,"a")
    file.write('epsilon' + ',' + 'Average Cost' + ',' + 'Cost variance' + ',' + 'Average Time' + '\n' )

def Jx_func(j,J,X):

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

def Jxx_func(j,J,X):

    #computing numberical hessians

    if j!=0 and j!=(X.size-1):
        Jxx = 0.5*((Jx_func(j+1,J,X) - Jx_func(j,J,X))/(X[j+1] - X[j]) +
               (Jx_func(j,J,X) - Jx_func(j-1,J,X))/(X[j] - X[j-1]))

    elif j==0:
        Jxx = (Jx_func(j+1,J,X) - Jx_func(j,J,X))/(X[j+1] - X[j])

    elif j==(X.size-1):
        Jxx = (Jx_func(j,J,X) - Jx_func(j-1,J,X))/(X[j] - X[j-1])

    return Jxx

def solve_DP(X, x0, n, N, epsilon = 0.0):

    dt = 1.0/N

    J = np.zeros((N+1,n)) #Value function
    J_temp = np.zeros(n)
    u = np.zeros((N,n)) #control

    #solving DP
    for i in range(N,-1,-1):

        if i == N:

            for j in range(X.size):
                J[i,j] = terminal_cost(X[j])

        else:

            for j in range(X.size):

                #print("i:", i, " j", j)
                Jx = Jx_func(j,J[i+1,:],X)
                u[i,j] = - Jx*g(X[j])/R #optimal control

                Jxx = Jxx_func(j,J[i+1,:],X)

                J[i,j] = current_cost(X[j],u[i,j],dt) + J[i+1,j] + Jx*g(X[j])*u[i,j]*dt +\
                            Jx*f(X[j])*dt + 0.5*Jxx*(epsilon**2)*dt

                Courant_num = abs(f(X[j])*dt/(X[1] - X[0]))
                #print('C num = ', Courant_num)

                if Courant_num >=1:
                    print('CFL condition violated!')


    #print "J:\n",J
    #print "u:\n",u

    return J,u

def monteCarloSims(X,x0,N,dt,epsi_range,iters=50):

    x_sol = np.zeros(N+1)
    x_sol[0] = x0
    x_sol[0] = X[find_nearest(X,x_sol[0])]

    for epsilon in epsi_range:

        J, u = solve_DP(X,x0,n,N,epsilon)

        cost_vector = np.zeros(iters)

        for iter in range(iters):
            U_opti = np.zeros(N)

            for i in range(N):

                idx = find_nearest(X,x_sol[i])
                #print idx
                try:
                    #print idx.size
                    U_opti[i] = u[i,idx[0]]

                except:

                    U_opti[i] = u[i,idx]

                x_sol[i+1] = model(x_sol[i],U_opti[i], dt, epsilon)


            #print "Solution:", x_sol
            #print "U:", U_opti
            cost = calculate_cost(x_sol,U_opti,dt)
            cost_vector[iter] = cost

        print "epsilon:",epsilon , "   Average Cost:",np.mean(cost_vector), "  Cost var:", np.var(cost_vector)

        if FILE_WRITE:
            file.write(str(epsilon) + ',' + str(np.mean(cost_vector)) + ',' + str(np.var(cost_vector)) + '\n')

def sampleTrial(X,x0,N,dt,epsilon = 0.0):

    J, u = solve_DP(X,x0,n,N,epsilon)

    x_sol = np.zeros(N+1)
    x_sol[0] = x0
    x_sol[0] = X[find_nearest(X,x_sol[0])]
    U_opti = np.zeros(N)

    for i in range(N):

        idx = find_nearest(X,x_sol[i])
        #print idx
        try:
            #print idx.size
            U_opti[i] = u[i,idx[0]]

        except:

            U_opti[i] = u[i,idx]

        x_temp = model(x_sol[i],U_opti[i], dt, epsilon)
        print('x_temp:',x_temp)
        x_sol[i+1] = X[find_nearest(X,x_temp)]

    print "Solution:"
    print "U:", U_opti
    print "X:", x_sol


    cost = calculate_cost(x_sol,U_opti,dt)
    print("cost:", cost)

    cost_to_go = np.zeros(N+1)
    for i in range(N+1):
        j = find_nearest(X,x_sol[i])
        cost_to_go[i] = J[i,j]

    #print "cost:",cost
    plot_func(J,X,N,x_sol,cost_to_go,dt)

if __name__=='__main__':

    x0 = 0.777
    n = 100 #state space discretization size
    N = 1000 #number of time steps
    dt = 1.0/N

    X = np.linspace(-1,1,n) #space discretization

    sampleTrial(X,x0,N,dt,epsilon = 0.0)
    #monteCarloSims(X,x0,N,dt, epsi_range = np.linspace(0.0,1.0,11),iters=50)
