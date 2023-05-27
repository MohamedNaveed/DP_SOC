# Dynamic Programming for a simple 1D case.
# Author: Mohamed Naveed Gul Mohamed
# email:mohdnaveed96@gmail.com
# Date: June 12th 2020
# Edited: Feb 21st 2021
from __future__ import division

import math as m
import numpy as np

import pylab
from mpl_toolkits import mplot3d #import Axes3D
from matplotlib import cm

params = {'axes.labelsize':12,
            'font.size':12,
            'legend.fontsize':14,
            'xtick.labelsize':10,
            'ytick.labelsize':10,
            'text.usetex':False,
            'figure.figsize':[5,5]}
pylab.rcParams.update(params)

R = 1.0
Q = 100.0
Qf = 500.0

u_max = 1#5 for 1D #1 for 1D cos

def current_cost(x, u, xT, dt):

    J = 0.5*(Q*((x - xT)**2) + R*(u**2))*dt

    return J

def terminal_cost(x, xT):

    J = 0.5*Qf*((x - xT)**2)

    return J

def f(x):

    #return(-x -2*(x**2) -0.5*(x**3))
    #return(-m.cos(x))
    return(x)

def g(x):

    return(1)

def model(x, u, dt, epsilon = 0.0):

    x_dot = f(x) + g(x)*u

    x_new = x + x_dot*dt + epsilon*u_max*m.sqrt(dt)*np.random.normal(0, 1, 1)

    return x_new

def find_u(x, x_next, dt):

    u = ((x_next - x)/dt - f(x))/g(x)

    return u

def plot_control(u_dp, u_lqr, x_dp, x_lqr, N, dt):
    
    pylab.figure(1)
    t = np.linspace(0,(N-1)*dt,N)
    pylab.plot(t,u_dp - u_lqr,'r.-')
    pylab.xlabel('Time')
    pylab.ylabel('difference in U')
    pylab.savefig('difference_u.png')
    
    pylab.figure(2)
    t = np.linspace(0,(N)*dt,N+1)
    pylab.plot(t,x_dp - x_lqr,'r.-')
    pylab.xlabel('Time')
    pylab.ylabel('difference in x')
    pylab.savefig('difference_x.png')

    


def plot_func(V,X,N,x_sol,cost,U,dt):

        fig = pylab.figure(1)
        ax = pylab.axes(projection='3d')
        #pylab.xlim(0,40)
        #pylab.ylim(0,12)
        #ax.yaxis.set_ticks(np.arange(0,6,1))
        #ax.yaxis.set_ticklabels(np.array([0,0.2,0.4,0.6,0.8,1]))
        #ax.xaxis.set_ticks(np.arange(0,5,1))
        #ax.xaxis.set_ticklabels(np.array([-1,-0.5,0,0.5,1]))
        #ax.zaxis.set_ticks(np.arange(0,4,1))
        #ax.xaxis.set_ticklabels(np.arange(0,40,5))
        pylab.grid()
        t = np.linspace(0,N*dt,N+1)
        x, y = np.meshgrid(X,t)

        # print "t:",t
        # print "x:",x
        # print "y:",y
        # print "V:",V
        t_steps =  np.linspace(0,N*dt,x_sol.size) #model time-steps

        ax.plot_wireframe(x, y, V)
        surf = ax.plot_surface(x, y, V,cmap=cm.RdYlGn_r) #plotting the cost-to-go surface
        ax.plot3D(x_sol,t_steps,cost,'r.-',markersize=5) #plotting the actual trajectory
        ax.set_ylabel('Time')
        ax.set_xlabel('X(0)')

        ax.set_zlabel('Cost-to-go')
        ax.set_zlim(0,500)
        fig.colorbar(surf,shrink=.3, aspect=10)
        #pylab.show()

        #plotting trajectory.
        pylab.figure(2)
        pylab.plot(t_steps,x_sol,'r.-')
        pylab.xlabel('Time')
        pylab.ylabel('X')
        #pylab.ylim(-0.2,1.1)

    
        pylab.figure(3)
        pylab.plot(t_steps[0:x_sol.size-1],U,'b.-')
        pylab.xlabel('Time')
        pylab.ylabel('U')
        pylab.show()

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def calculate_cost(x, u, xT, dt):

    N = u.size
    cost = 0
    for i in range(N):
        cost += current_cost(x[i],u[i],xT,dt)
        #print("cost: ", cost)
        #print("i:", i, " Cost:", round(cost,2), " x:",round(x[i],2), " u:", round(u[i],2))
    cost += terminal_cost(x[N],xT)
    #print("i:", N, " Cost:", round(cost,2), " x:", round(x[N],2))

    return cost

def Exp_J(x, u, J, X, dt, epsilon):

    #calculating expectation of J using monte carlo

    N_samples = 1 # number of monte carlo samples.

    cost = 0

    for i in range(N_samples):

        x_new = model(x, u, dt, epsilon) # next state
        x_newIdx = find_nearest(X, x_new) # finding closest point in the discretization

        cost += J[x_newIdx] # adding all the costs

        #print('x = ', x, ' x_new = ', x_new, " J = ", J[x_newIdx])

    #print('Cost = ',cost/N_samples)
    return cost/N_samples # average - expectation
