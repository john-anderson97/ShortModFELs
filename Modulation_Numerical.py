##FEL Motion Numerical (Runge Kutta 4th Order Method)
#John Anderson

#importation of required packages used in program
import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.integrate import solve_ivp
import matplotlib.gridspec as gridspec
import matplotlib
import sys
#formats the font for all matplotlib plots
font = {'family' : 'verdana', 'size'   : 14}
matplotlib.rc('font', **font)

start_time=time.time()  #create an initial time for when the code was run

#mass of electron (kg), charge of electron (C), magnetic field strength (T), speed of light (m/s),velocity of travelling particle (m/s)
m = 9.11e-31; q = 1.6e-19; B = .25;c = 2.99792458e8; v=0.9999998*c; alpha = 0.1
beta = v/c #speed of the electron as a factor of the speed of light
lorentz_factor=1/(np.sqrt(1-beta**2))
k = (np.pi * 2) / 2e-2  # wave number defined as 2pi divided by the wavelength (1/m)
au =(q * B) / (k * m * c) #undulator parameter
n =2

def f1(t,y):    #x direction coupled equations
    # [x' = u, x'' = u']
    return [y[1], (q*B*beta*c)/(lorentz_factor*m)*(np.cos(k*beta*c*t)*(1+alpha*np.sin(n*k*beta*c*t)))]

#Short Modulation greater than or equal to 1.
if n >=0:
    #Equal Modulation Case
    if n == 1:
        print(n)
        def f2(t, y):  # z direction equations
            # [z' = u]
            r = k * beta * c * t
            return [(beta * c) + ((c * au ** 2) / (4 * lorentz_factor ** 2)) * (
                        np.cos(2 * r) + alpha * ((np.sin(r) * np.cos(2 * r)) - ((alpha * np.cos(4 * r)) / 16)))]
        # Initial conditions
        x0 = -au / (lorentz_factor * k * beta)
        vx0 = -(q * B) / (lorentz_factor * m * k) * (alpha / 4)
        z0 = ((alpha * au ** 2) / (12 * (lorentz_factor ** 2) * beta * k))
    #Half Modulation Case
    elif n==2:
        def f2(t, y):  # z direction equations
            # [z' = u]
            r = k * beta * c * t

            d1 = (c * au ** 2) / (4 * lorentz_factor ** 2)
            d2 = np.cos(2 * r)
            d3 = alpha * np.sin(2 * r)
            d4 = (alpha / 3) * np.sin(4 * r)
            d5 = (alpha / 3) * np.sin(2 * r)
            d6 = alpha * np.cos(r)
            d7 = (alpha / 3) * np.cos(3 * r)

            return [beta * c + d1 * (d2 + d3 + d4 - d5 - d6 - d7)]
            # Initial conditions
        x0 = -au / (lorentz_factor * k * beta)
        vx0 = -(2 * alpha * q * B) / (3 * lorentz_factor * m * k)
        z0 = -(5 * alpha * au ** 2) / (96 * k * beta * lorentz_factor ** 2)
    #General Modulation Case
    else:
        print(n)
        def f2(t, y):  # z direction equations
            # [z' = u]
            r = k * beta * c * t

            d1 = (c * au ** 2) / (4 * lorentz_factor ** 2)
            d2 = np.cos(2 * r)
            d3 = (alpha / (n + 1)) * (np.sin((n + 2) * r) + np.sin(n * r))
            d4 = (alpha / (n - 1)) * (-np.sin((n - 2) * r) + np.sin(n * r))
            d5 = ((alpha / (2 * (n + 1))) ** 2) * np.cos(2 * (n + 1) * r)
            d6 = ((alpha ** 2) / (2 * (n + 1) * (n - 1))) * (np.cos(2 * r) + np.cos(2 * n * r))
            d7 = ((alpha / (2 * (n - 1))) ** 2) * np.cos(2 * (n - 1) * r)
            return [beta * c + d1 * (d2 + d3 + d4 - d5 - d6 - d7)]
        x0 = -au / (lorentz_factor * k * beta)
        vx0 = (-(alpha * q * B) / (2 * lorentz_factor * m * k)) * ((1 / (n + 1)) + (1 / (n - 1)))
        z0 = ((alpha * au ** 2) / (2 * k * beta * lorentz_factor ** 2)) * ((4 - n) / ((n + 2) * (n - 1) * (n - 2)))
else:
    print('n must have a value greater than or equal to 1')
    sys.exit()

#RK45 method plus plotting
def main():
    tTime = 1e-9  # total time
    stepsize = tTime / 5000  # number of steps

    ##solutions to equations using RK45 method with initial condidtion, timeframe and step size.
    x_sol = solve_ivp(f1, t_span=[0, tTime], y0=[x0, vx0], method='RK45', max_step=stepsize, atol=1, rtol=1)
    z_sol = solve_ivp(f2, t_span=[0, tTime], y0=[z0], method='RK45', max_step=stepsize, atol=1, rtol=1)


    t1 = z_sol.t  # data for time, same for x and z
    x = x_sol.y[0]  # data for x
    vx = x_sol.y[1]  # data for vx
    z = z_sol.y[0]  # data for z
    vz = np.squeeze(f2(t1, 0))  # data for vz
    plt.plot(t1, vz)
    print(max(vz))
    print(max(vx))
    plt.show()
    plt.close()
    scale_x = x / 9.385366059602438e-07  # scaled x
    scale_time1 = t1 / 1e-9  # scaled time
    scale_z = 3.459093113508857e-11

    gs = gridspec.GridSpec(4, 2)
    fig = plt.figure(figsize=(10, 5))
    plt.suptitle(r'$\lambda_m = \lambda_u /$'+str(n))
    # plot x(t) against time
    plt.subplot(gs[0, :], xlabel='t', ylabel='x(t)')
    plt.plot(scale_time1, scale_x, 'red')
    # plot z(t) minus beta*c*t against time
    plt.subplot(gs[1, :], xlabel='t', ylabel=r'$z(t)-\beta ct$')
    plt.plot(scale_time1, ((z - (beta * c * t1)) / scale_z), 'blue')
    # plot z(t) minus beta*c*t against x(t)
    plt.subplot(gs[2:, :], xlabel=r'$z(t)-\beta ct$', ylabel='x(t)')
    plt.plot((z - (beta * c * t1)) / scale_z, scale_x, 'purple')
    matplotlib.pyplot.subplots_adjust(wspace=.75, hspace=.75)
    plt.show()

    gs = gridspec.GridSpec(3, 2)
    plt.figure(figsize=(5, 5))
    plt.suptitle(r'$\lambda_m = \lambda_u /$' + str(n))
    # calculation of total speed
    v = ((vx ** 2) + (vz ** 2)) ** .5

    # plot vx(t) against scaled time
    plt.subplot(gs[0, :], xlabel='t', ylabel=r'$\beta_x$(t)')
    plt.plot(scale_time1, vx/88393.76523535939 , 'red')

    # plot vz(t) against scaled time
    plt.subplot(gs[1, :], xlabel='t', ylabel=r'$\beta_z$(t)')
    plt.plot(scale_time1, vz/299792404.5572332 , 'blue')
    c2 = (299792404.5572332**2+88393.76523535939**2)**0.5
    # plot v(t) against scaled time
    plt.subplot(gs[2, :], xlabel='t', ylabel=r'$\beta$(t)')
    plt.plot(scale_time1, v/c , 'purple')

    matplotlib.pyplot.subplots_adjust(wspace=.75, hspace=.75)
    plt.show()
    plt.close()
    print(str(round((time.time() - start_time), 4)) + 's')

main()






