##FEL Motion Numerical (Runge Kutta 4th Order Method)
#John Anderson

#importation of required packages used in program

import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.integrate import solve_ivp
import matplotlib.gridspec as gridspec
import matplotlib

#formats the font for all maplotlib plots
font = {'family' : 'verdana', 'size'   : 14}
matplotlib.rc('font', **font)

start_time=time.time()  #create an initial time for when the code was run

m = 9.11e-31 #mass of electron (kg)
q = 1.6e-19 #charge of electron (C)
B = .25 #magnetic field strength (T)

c = 2.99792458e8 #speed of light (m/s)
v=0.9999998*c #velocity of travelling particle with speed of % times the speed of light (m/s)
beta = v/c #speed of the electron as a factor of the speed of light
lorentz_factor=1/(np.sqrt(1-beta**2))
k = (np.pi * 2) / 2e-2  # wave number defined as 2pi divided by the wavelength (1/m)
au =(q*B)/(m*c*k)
alpha = 0.1

def f1(t,y):    #x direction coupled equations
    # x' = u
    x1 = (q*B*beta*c)/(lorentz_factor*m)
    x2 = np.cos(k*beta*c*t)
    x3 = (1+alpha*np.sin(k*beta*c*t))
    return [y[1], x1*x2*x3]

def f2(t,y):    #z direction coupled equations
    x1 = (au**2)*(np.cos(2*k*beta*c*t)/2)
    x2 = (((alpha*au)/2)**2)*((np.cos(4*k*beta*c*t))/2)
    x3 = (c/(2*(lorentz_factor**2)))
    return [beta*c +(x3*x1)-(x3*x2)]

# Initial conditions
x0 = -au / (lorentz_factor * k * beta)
vx0 = ((q * B) / (lorentz_factor * m * k)) * (-alpha / 4)
z0 = 0

tTime = 1e-9  # total time
stepsize = tTime / 10000  # number of steps

##solutions to equations using RK45 method with initial condidtion, timeframe and step size.
x_sol = solve_ivp(f1, [0, tTime], [x0, vx0], method='RK45', max_step=stepsize)
z_sol = solve_ivp(f2, [0, tTime], [z0], method='RK45', max_step=stepsize)

t1 = z_sol.t  # data for time, same for x and z

x = x_sol.y[0]  # data for x
vx = x_sol.y[1]  # data for vx

z = z_sol.y[0]  # data for z
vz = np.squeeze(f2(t1, 0))  # data for vz

scale_x = 9.385366059602438e-07
scale_z = 3.459093113508857e-11
scale_time = 1e-09
gs = gridspec.GridSpec(4, 2)
fig = plt.figure(figsize=(10, 5))

plt.subplot(gs[0, :], xlabel='t', ylabel='x(t)')
plt.plot(t1/scale_time, x/scale_x, 'red')

plt.subplot(gs[1, :], xlabel='t', ylabel='z(t)')
plt.plot(t1/scale_time, ((z - (beta * c * t1))/scale_z), 'blue')

plt.subplot(gs[2:, :], xlabel=r'$z(t)-\beta ct$', ylabel='x(t)')
plt.plot(((z - (beta * c * t1))/scale_z), x/scale_x, 'purple')
matplotlib.pyplot.subplots_adjust(wspace=.75, hspace=.75)

fig = plt.figure(figsize=(10, 5))
plt.plot(t1, vx, 'red')
plt.plot(t1, vz, 'blue')
v = (vx ** 2) + (vz ** 2)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
# plt.ylim(0,c+1e8)
beta_2 = (v ** .5) / c
lorentz = ((1 / (1 - (beta_2 ** 2)) ** .5))
print(np.mean(lorentz))


# prints the difference of the start time to finish time, formatted to show '-time (s) -'
print("-%.2fs-" % (time.time() - start_time))
# shows the plots in the window
plt.show()
## Time to run code: 0.58s


