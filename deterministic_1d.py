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
    n = 10 #state space discretization size
    X = np.linspace(-1,1,n)

    N = 100 #number of time steps
    dt = 1.0/N
    J = np.zeros((N+1,n)) #Value function
    J_temp = np.zeros(n)
    u = np.zeros((N,n)) #control

    x_sol = np.zeros(N+1)

    for i in range(N,-1,-1):

        if i == N:

            for j in range(X.size):
                J[i,j] = terminal_cost(X[j])

        else:

            for j in range(X.size):
                for k in range(X.size):
                    J_temp[k] = current_cost(X[j],find_u(X[j],X[k],dt),dt) + J[i+1,k]

                x_next_idx = np.argmin(J_temp)
                u[i,j] = find_u(X[j],X[x_next_idx],dt)
                J[i,j] = current_cost(X[j],u[i,j],dt) + J[i+1,x_next_idx]


    print "X:\n",X
    print "J:\n",J
    print "u:\n",u

    x_sol[0] = x0
    x_sol[0] = X[find_nearest(X,x_sol[0])]
    U_opti = np.zeros(N)

    for i in range(N):

        idx = find_nearest(X,x_sol[i])

        try:
            #print idx.size
            U_opti[i] = u[i,idx[0]]

        except:

            U_opti[i] = u[i,idx]

        x_sol[i+1] = model(x_sol[i],U_opti[i], dt)

    print "Solution:", x_sol
    print "U:", U_opti
    cost = calculate_cost(x_sol,U_opti,dt)
    print("cost:", cost)

    cost_to_go = np.zeros(N+1)
    for i in range(N+1):
            j = find_nearest(X,x_sol[i])
            cost_to_go[i] = J[i,j]

    print "cost:",cost_to_go
    plot_func(J,X,N,x_sol,cost_to_go,dt)
