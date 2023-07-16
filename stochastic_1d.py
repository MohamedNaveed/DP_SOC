# Dynamic Programming for a simple 1D case.
# Author: Mohamed Naveed Gul Mohamed
# email:mohdnaveed96@gmail.com
# Date: Feb 21st 2021
from __future__ import division

from dp_funcs import *
import math as m
import numpy as np
from scipy.io import savemat
import time

FILE_WRITE = False# False | True

if FILE_WRITE:
    filename = "stochastic_hjb_1dcos_T300000_X200_processNoise_e0.csv"
    file = open(filename,"a")
    file.write('epsilon' + ',' + 'Average Cost' + ',' + 'Cost variance' + ',' + 'Time taken' + '\n' )

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

def solve_DP_Bellman_v2(X, x0, xT, n, N, dt, epsilon = 0.0):


    # use optimal control to transition between states. 

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
                
                #Calculating the control $u$ to take X(j) to X(k)
                Jx = Jx_func(j,J[t+1,:],X)
                u[t,j] = - Jx*g(X[j])/R #optimal control
                
                x_new = model(X[j], u[t,j], dt, epsilon) # next state
                x_newIdx = find_nearest(X, x_new) # finding closest point in the discretization

                   
                J[t,j] = current_cost(X[j], u[t,j], xT, dt) + J[t+1, x_newIdx]

                #print("t: ", t, " j: ", j, " J:", J[t,j], " u : ", u[t,j], "Jx :", Jx)


    return J,u


def solve_DP_Bellman(X, x0, xT, n, N, dt, epsilon = 0.0):

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

        print("t = ", t)

    #print "J:\n",J
    #print "u:\n",u

    return J,u

def solve_DP_HJB(X, x0, xT, n, N, dt, epsilon = 0.0):


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

def solve_LQR(N_steps, del_t, epsilon):
    
    K = np.zeros(N_steps) # feedback gain

    P = np.zeros(N_steps+1) # cost_to_go state
    q = np.zeros(N_steps+1) # cost_to_go linear term due to noise

    Ac = 1;#linear system A continuous
    Bc = 1;# linear system B continuous
    
    Ad = Ac + 1*del_t #linear system A discrete
    Bd = Bc*del_t # linear system B discrete
    

    setting = "continuous"; #"discrete","continuous"

    if setting == "discrete":
        
        #continuous time lqr
        for i in range(N_steps,-1,-1):
        
            if i == N_steps:
                P = Qf;
            
            else:
                P_dot = -(Ac*P + P*Ac - (P*Bc*Bc*P/R) + Q)
                
                P = P + P_dot*del_t
                K[i] = -(1.0/R)*Bc*P

                #print("K:", K[i])

        return K, P, Ac, Bc
    
        
        
    else:
        #discrete time lqr
        print("executing discrete")
        for i in range(N_steps,-1,-1):
            
            if i == N_steps:
                
                P[i] = Qf;
                q[i] = 0;
                
            else:
                
                K[i] = - (1.0/(R*del_t + Bd*P[i+1]*Bd))*Bd*P[i+1]*Ad
                #print("K:", K[i])
                P[i] = Q*del_t + Ad*P[i+1]*Ad + Ad*P[i+1]*Bd*K[i]
                
                q[i] = q[i+1] + 0.5*P[i+1]*(epsilon**2)*del_t
            
        
        return K,P,q,Ad,Bd

    
    
def monteCarloSims(X, x0, xT, N, dt, del_t, FIX_DP_EPSILON, epsi_range, iters=50):

    N_steps = int(N*dt / del_t) # number of steps the system propagates

    factor_del_t = int(del_t / dt)

    x_sol = np.zeros(N_steps+1)
    x_sol[0] = x0
    
    x_sol_lqr = np.zeros(N_steps+1)
    x_sol_lqr[0] = x0

    #J, u = solve_DP_Bellman(X, x0, xT, n, N, 0.0)
    #J, u = solve_DP_HJB(X, x0, xT, n, N, 0.0)
    if FIX_DP_EPSILON:
        J, u = solve_DP_HJB(X, x0, xT, n, N, dt, 0.0)
    
    for epsilon in epsi_range:
        
        print("epsilon:", epsilon)
        start_exec = time.time()
        
        if not FIX_DP_EPSILON:
            J, u = solve_DP_HJB(X, x0, xT, n, N, dt, epsilon)
        #J, u = solve_DP_Bellman(X, x0, xT, n, N, epsilon)
        #LQR solution
        
        #K, P, q, _, _ = solve_LQR(N, dt, epsilon)
    
        cost_vector = np.zeros(iters)
        #cost_vector_lqr = np.zeros(iters)

        for iter in range(iters):
            print("iter:", iter)
            U_opti = np.zeros(N_steps)
            #U_opti_lqr = np.zeros(N_steps)

            for i in range(N_steps):
                
                #DP
                idx = find_nearest(X,x_sol[i])

                U_opti[i] = u[i*factor_del_t,idx]

                x_sol[i+1] = model(x_sol[i],U_opti[i], del_t, epsilon)
                
                #LQR
                #U_opti_lqr[i] = -K[i*factor_del_t]*(xT - x_sol_lqr[i])

                #x_sol_lqr[i+1] = model(x_sol_lqr[i],U_opti_lqr[i], del_t, epsilon)

            #print "Solution:", x_sol
            #print "U:", U_opti
            cost = calculate_cost(x_sol, U_opti, xT, del_t)
            cost_vector[iter] = cost
            
            #cost_lqr = calculate_cost(x_sol_lqr, U_opti_lqr, xT, del_t)
            #cost_vector_lqr[iter] = cost_lqr

        end_exec = time.time()
        time_taken = (end_exec - start_exec)/iters
        print("epsilon:", epsilon , "   Average Cost:",np.mean(cost_vector), "  Cost var:", np.var(cost_vector))
        #print("LQR: epsilon:", epsilon , "   Average Cost:",np.mean(cost_vector_lqr), "  Cost var:", np.var(cost_vector_lqr))
        
        if FILE_WRITE:
            file.write(str(epsilon) + ',' + str(np.mean(cost_vector)) + ',' + str(np.var(cost_vector)) + ',' + str(time_taken) + '\n')
            #file.write('LQR' + ',' + str(epsilon) + ',' + str(np.mean(cost_vector_lqr)) + ',' + str(np.var(cost_vector_lqr)) + '\n')

def sampleTrial(X, x0, xT, N, dt, del_t, epsilon = 0.0):
    
    N_steps = int(N*dt / del_t) # number of steps the system pro

    factor_del_t = int(del_t / dt)

    #print('N_steps = ', N_steps, ' Factor = ', factor_del_t)
    

    #J, u = solve_DP_Bellman_v2(X, x0, xT, n, N, dt, epsilon)

    J, u = solve_DP_HJB(X, x0, xT, n, N, dt, epsilon)
    
    #J, u = solve_DP_Bellman(X, x0, xT, n, N, dt, epsilon)

    print('DP solved.')
    
    x_sol = np.zeros(N_steps+1)
    x_sol[0] = x0
    U_opti = np.zeros(N_steps)

    print('Execution epsilon:', epsilon)

    for i in range(N_steps):

        idx = find_nearest(X,x_sol[i])

        U_opti[i] = u[i*factor_del_t, idx]

        x_sol[i+1] = model(x_sol[i],U_opti[i], del_t, epsilon)

    print("Solution:")
    #for t, x_temp, u_temp in zip(range(N_steps),x_sol, U_opti):
    #    print('t =', t, '   X:', x_temp, '   U:',u_temp)
    print('t =', N_steps, '   X:', x_sol[N_steps])

    cost = calculate_cost(x_sol,U_opti,xT,del_t)
    print("cost:", cost)

    cost_to_go = np.zeros(N_steps+1)
    for i in range(N_steps+1):
        j = find_nearest(X,x_sol[i])
        cost_to_go[i] = J[i*factor_del_t,j]
    
    print("cost-to-go DP:", cost_to_go[0])
    #plot_func(J,X,N,x_sol,cost_to_go,U_opti,dt)
    
    
    #LQR solution
    '''
    K, P, q, _, _ = solve_LQR(N, dt, epsilon)
    x_sol_lqr = np.zeros(N_steps+1)
    x_sol_lqr[0] = x0
    U_opti_lqr = np.zeros(N_steps)

    for i in range(N_steps):

        U_opti_lqr[i] = -K[i*factor_del_t]*(xT - x_sol_lqr[i])

        x_sol_lqr[i+1] = model(x_sol_lqr[i],U_opti_lqr[i], del_t, epsilon)
    
    print("Solution LQR:")
    #for t, x_temp, u_temp in zip(range(N_steps),x_sol_lqr, U_opti_lqr):
        #print('t =', t, '   X:', x_temp, '   U:',u_temp)
    print('t =', N_steps, '   X:', x_sol_lqr[N_steps])

    cost_lqr = calculate_cost(x_sol_lqr,U_opti_lqr,xT,del_t)
    print("cost:", cost_lqr)
    
    #cost_to_go_lqr = 0.5*P[0]*((xT - x0)**2)  
    #print("cost-to-go lqr without q:", cost_to_go_lqr)
    
    cost_to_go_lqr = 0.5*P[0]*((xT - x0)**2) + q[0]
    print("cost-to-go lqr with q:", cost_to_go_lqr)
    
    diff_u = U_opti - U_opti_lqr
    diff_x = x_sol - x_sol_lqr
    #plot_control(U_opti, U_opti_lqr, x_sol, x_sol_lqr, N_steps, del_t)
    var_dic = {"u_dp": U_opti, "u_lqr": U_opti_lqr, 
                "x_dp": x_sol, "x_lqr": x_sol_lqr, "N": N, "dt": dt,
                "J": J, "u_global": u, "X": X, "P": P, "K": K, "q": q}
    savemat("dp_hjb_lqr_n_200_N_300000_epsi1000.mat",var_dic)
    '''
    var_dic = {"u_dp": U_opti, "x_dp": x_sol, "N": N, "dt": dt,
                "J": J, "u_global": u, "X": X}
    savemat("dp_hjb_1dpoly_n_200_N_300000_epsi0.mat",var_dic)
    
def check_LQR(N, dt, x0, xT):

    #LQR solution
    print("solving LQR ....")
    K, P,Ad,Bd = solve_LQR(N, dt)
    x_sol_lqr = np.zeros(N+1)
    x_sol_lqr[0] = x0
    U_opti_lqr = np.zeros(N)

    for i in range(N):

        U_opti_lqr[i] = -K[i]*(xT - x_sol_lqr[i])
        #x_sol_lqr[i+1] = Ad*x_sol_lqr[i] + Bd*U_opti_lqr[i]
        x_sol_lqr[i+1] = model(x_sol_lqr[i],U_opti_lqr[i], dt)

    print("Solution LQR:")
    for t, x_temp, u_temp in zip(range(N),x_sol_lqr, U_opti_lqr):
        print('t =', t, '   X:', x_temp, '   U:',u_temp)
    print('t =', N, '   x:', x_sol_lqr[N])

    cost_lqr = calculate_cost(x_sol_lqr,U_opti_lqr,xT,del_t)
    print("cost:", cost_lqr)
    
    cost_to_go_lqr = 0.5*P[0]*((xT - x0)**2)
    print("cost-to-go lqr:", cost_to_go_lqr)

    K, P = solve_LQR(N, dt)
    x_sol_lqr = np.zeros(N_steps+1)
    x_sol_lqr[0] = x0
    U_opti_lqr = np.zeros(N_steps)

    for i in range(N_steps):

        U_opti_lqr[i] = -K[i*factor_del_t]*(xT - x_sol_lqr[i])

        x_sol_lqr[i+1] = model(x_sol_lqr[i],U_opti_lqr[i], del_t, epsilon)
    
    print("Solution LQR:")
    #for t, x_temp, u_temp in zip(range(N_steps),x_sol_lqr, U_opti_lqr):
        #print('t =', t, '   X:', x_temp, '   U:',u_temp)
    print('t =', N_steps, '   X:', x_sol_lqr[N_steps])

    cost_lqr = calculate_cost(x_sol_lqr,U_opti_lqr,xT,del_t)
    print("cost:", cost_lqr)
    
    cost_to_go_lqr = 0.5*P*((xT - x0)**2)
    print("cost-to-go lqr:", cost_to_go_lqr)
    
    plot_control(U_opti[i], U_opti_lqr, x_sol, x_sol_lqr, N, dt)
    
if __name__=='__main__':

    x0 = 1.0 #initial state.

    xT = 0.0 #final state.

    n = 200 #state space discretization size
    N = 300000 #number of time steps
    dt = 1.0/N #dt - time step for DP solution.

    del_t = dt #0.01 #del_t - time step for model update.

    X = np.linspace(-2,2,n) #space discretization

    #np.random.seed(2)

    sampleTrial(X, x0, xT, N, dt, del_t, epsilon = 0.0)
    #check_LQR(N, dt, x0, xT)
    FIX_DP_EPSILON = True
    epsi_range = np.array([0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,4.0,6.0,8.0,10.0])
    #monteCarloSims(X, x0, xT, N, dt, del_t, FIX_DP_EPSILON, epsi_range,iters=100)
