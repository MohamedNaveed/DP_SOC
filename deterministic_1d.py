# Dynamic Programming for a simple 1D case.
# Author: Mohamed Naveed Gul Mohamed
# email:mohdnaveed96@gmail.com
# Date: June 12th 2020
from __future__ import division

from dp_funcs import *
import math as m
import numpy as np

if __name__=='__main__':

    x0 = 0.8
    n = 1000 #state space discretization size
    X = np.linspace(-1,1,n)

    N = 100 #number of time steps
    dt = 1.0/N
    J = np.zeros((N+1,n)) #Value function
    J_temp = np.zeros(n)
    u = np.zeros((N,n)) #control

    x_sol = np.zeros(N+1)
    print "X:\n",X

    for t in range(N,-1,-1):

        if t == N:

            for j in range(X.size):
                J[t,j] = terminal_cost(X[j])

            #print('t = ', t, ' J = ', J[t,:])

        else:

            for j in range(X.size):
                for k in range(X.size):

                    u_temp = find_u(X[j],X[k],dt)
                    J_temp[k] = current_cost(X[j], u_temp, dt) + J[t+1,k]

                    #print('k = ', k, 'u =', u_temp, ' inc cost = ', current_cost(X[j],find_u(X[j],X[k],dt),dt), ' J(t+1,k) = ', J[t+1,k])





                x_next_idx = np.argmin(J_temp)


                #print('j = ', j, '   J_temp', J_temp)
                #print('current idx = ', j,'  x_next_idx = ', x_next_idx)

                u[t,j] = find_u(X[j],X[x_next_idx],dt)
                J[t,j] = current_cost(X[j],u[t,j],dt) + J[t+1,x_next_idx]

            #print('t = ', t, ' J = ', J[t,:])


    #print "u:\n",u

    x_sol[0] = x0
    U_opti = np.zeros(N)

    for i in range(N):

        idx = find_nearest(X,x_sol[i])

        try:
            #print idx.size
            U_opti[i] = u[i,idx[0]]

        except:

            U_opti[i] = u[i,idx]

        x_sol[i+1] = model(x_sol[i], U_opti[i], dt)


    print "Solution:", x_sol
    print "U:", U_opti
    cost = calculate_cost(x_sol,U_opti,dt)
    print("cost:", cost)

    cost_to_go = np.zeros(N+1)
    for i in range(N+1):
            j = find_nearest(X,x_sol[i])
            cost_to_go[i] = J[i,j]

    #print "cost to go:",cost_to_go
    plot_func(J,X,N,x_sol,cost_to_go,dt)
