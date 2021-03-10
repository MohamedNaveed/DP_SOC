# Dynamic Programming for a simple 1D case.
# Author: Mohamed Naveed Gul Mohamed
# email:mohdnaveed96@gmail.com
# Date: June 12th 2020
# Edited: Feb 21st 2021

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
            'text.usetex':True,
            'text.latex.unicode':True,
            'figure.figsize':[4.5,4.5]}
pylab.rcParams.update(params)

R = 1.0
Q = 1.0
Qf = 10.0

def current_cost(x,u,dt):

    J = 0.5*(Q*(x**2) + R*(u**2))*dt

    return J

def terminal_cost(x):

    J = 0.5*Qf*(x**2)

    return J

def f(x):

    return(-x -2*(x**2) -0.5*(x**3))

def g(x):

    return(1)

def model(x, u, dt, epsilon = 0.0):

    x_dot = f(x) + g(x)*u

    x_new = x + x_dot*dt + epsilon*m.sqrt(dt)*np.random.normal(0, 1, 1)

    return x_new

def find_u(x, x_next, dt):

    u = (x_next - x)/dt - (-x -2*(x**2) -0.5*(x**3))

    return u

def plot_func(V,X,N,x_sol,cost,dt):

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

        ax.plot_wireframe(x, y, V)
        surf = ax.plot_surface(x, y, V,cmap=cm.RdYlGn_r) #plotting the cost-to-go surface
        ax.plot3D(x_sol,t,cost,'r.-',markersize=5) #plotting the actual trajectory
        ax.set_ylabel('Time')
        ax.set_xlabel('X(0)')

        ax.set_zlabel('Cost-to-go')
        fig.colorbar(surf,shrink=.3, aspect=10)
        #pylab.show()

        #plotting trajectory.
        pylab.figure(2)

        pylab.plot(t,x_sol,'r.-')
        pylab.xlabel('Time')
        pylab.ylabel('X')
        pylab.ylim(-0.,1.1)
        pylab.show()

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def calculate_cost(x,u,dt):

    N = u.size
    cost = 0
    for i in range(N):
        cost += current_cost(x[i],u[i],dt)
        #print("i:", i, " Cost:", cost, " x:",x[i], " u:", u[i])
    cost += terminal_cost(x[N])
    #print("i:", N, " Cost:", cost, " x:",x[N])
    return cost
