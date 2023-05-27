from __future__ import division

import math as m
import numpy as np
from scipy.io import loadmat
import pylab

params = {'axes.labelsize':12,
            'font.size':12,
            'legend.fontsize':14,
            'xtick.labelsize':10,
            'ytick.labelsize':10,
            'text.usetex':False,
            'figure.figsize':[5,5]}
pylab.rcParams.update(params)

def plot_data(u_dp, u_lqr, x_dp, x_lqr, N, dt):
    
    t_u = np.linspace(0,(N-1)*dt,N).reshape(1,N)
    t_x = np.linspace(0,(N)*dt,N+1).reshape(1,N+1)

    
    pylab.figure(1)
    
    pylab.plot(t_u,u_dp - u_lqr,'r--')
    pylab.xlabel('Time')
    pylab.ylabel('difference in U')
    pylab.savefig('dp_lqr/difference_u.png')
    
    pylab.figure(2)

    pylab.plot(t_x,x_dp - x_lqr,'r--')
    pylab.xlabel('Time')
    pylab.ylabel('difference in x')
    pylab.savefig('dp_lqr/difference_x.png')

    pylab.figure(3)
    
    pylab.plot(t_u,u_dp,'r.-')
    pylab.xlabel('Time')
    pylab.ylabel('U DP')
    pylab.savefig('dp_lqr/u_dp.png')

    pylab.figure(4)
    
    pylab.plot(t_u,u_lqr,'r.-')
    pylab.xlabel('Time')
    pylab.ylabel('U LQR')
    pylab.savefig('dp_lqr/u_lqr.png')

    pylab.figure(5)
    
    pylab.plot(t_x,x_dp,'r.-')
    pylab.xlabel('Time')
    pylab.ylabel('x DP')
    pylab.savefig('dp_lqr/x_dp.png')

    pylab.figure(6)
    
    pylab.plot(t_x,x_lqr,'r.-')
    pylab.xlabel('Time')
    pylab.ylabel('x LQR')
    pylab.savefig('dp_lqr/x_lqr.png')



if __name__=='__main__':

    data_dic = loadmat("dp_lqr.mat")

    u_dp = data_dic["u_dp"]
    u_lqr = data_dic["u_lqr"]
    x_dp = data_dic["x_dp"]
    x_lqr = data_dic["x_lqr"]

    N = np.size(u_dp)
    dt = 1.0/N

    plot_data(u_dp, u_lqr, x_dp, x_lqr, N, dt)