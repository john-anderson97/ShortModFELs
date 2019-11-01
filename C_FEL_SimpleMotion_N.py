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

#mass of electron (kg), charge of electron (C), magnetic field strength (T), speed of light (m/s),velocity of travelling particle (m/s)
m = 9.11e-31; q = 1.6e-19; B = .25;c = 2.99792458e8; v=0.9999998*c
beta = v/c #speed of the electron as a factor of the speed of light
lorentz_factor=1/(np.sqrt(1-beta**2))
k = (np.pi * 2) / 2e-2  # wave number defined as 2pi divided by the wavelength (1/m)
au =(q * B) / (k * m * c) #undulator parameter


def f1(t,y):    #x direction coupled equations
    # [x' = u, x'' = u']
    return [y[1], (q*B*beta*c)/(lorentz_factor*m)*np.cos(k*beta*c*t)]

def f2(t,y):    #z direction equations
    # [z' = u]
    return [beta*c +((c/(2*(lorentz_factor**2)))*(au**2)*(np.cos(2*k*beta*c*t)/2))]

# Initial conditions
x0 = -au / (lorentz_factor * k * beta)
vx0 = 0
z0 = 0

tTime = 1e-9 # total time
stepsize = tTime / 10000 # number of steps

##solutions to equations using RK45 method with initial condidtion, timeframe and step size.
x_sol = solve_ivp(f1, [0, tTime], [x0, vx0], method='RK45', max_step=stepsize,atol = 1, rtol = 1)
z_sol = solve_ivp(f2, [0, tTime], [z0], method='RK45', max_step=stepsize, atol = 1, rtol = 1)

t1 = z_sol.t  # data for time, same for x and z

x = x_sol.y[0]  # data for x
vx = x_sol.y[1]  # data for vx

z = z_sol.y[0]  # data for z
vz = np.squeeze(f2(t1, 0))  # data for vz

scale_x = x/max(x) #scaled x
scale_z = 3.459093113508857e-11 #scaled z
scale_time =t1/max(t1) #scaled time

gs = gridspec.GridSpec(4, 2)
fig = plt.figure(figsize=(10, 5))
#plot x(t) against time
plt.subplot(gs[0, :], xlabel='t', ylabel='x(t)')
plt.plot(scale_time, scale_x, 'red')
#plot z(t) minus beta*c*t against time
plt.subplot(gs[1, :], xlabel='t', ylabel=r'$z(t)-\beta ct$')
plt.plot(scale_time, ((z - (beta * c * t1))/scale_z), 'blue')
#plot z(t) minus beta*c*t against x(t)
plt.subplot(gs[2:, :], xlabel=r'$z(t)-\beta ct$', ylabel='x(t)')
plt.plot(((z - (beta * c * t1))/scale_z), scale_x, 'purple')
matplotlib.pyplot.subplots_adjust(wspace=.75, hspace=.75)

plt.figure(figsize=(5, 5))
#calculation of total speed
v = ((vx ** 2) + (vz ** 2))**.5

# plot vx(t) against scaled time
plt.subplot(gs[0, :], xlabel='t', ylabel=r'$\beta_x$(t)')
plt.plot(scale_time, vx/c, 'red')

# plot vz(t) against scaled time
plt.subplot(gs[1, :], xlabel='t', ylabel=r'$\beta_z$(t)')
plt.plot(scale_time, vz/c, 'blue')

# plot v(t) against scaled time
plt.subplot(gs[2, :], xlabel='t', ylabel=r'$\beta$(t)')
plt.plot(scale_time, v/c, 'purple')

plt.subplot(gs[3, :], xlabel='t')
#plot all velocities against scale time
plt.plot(scale_time, vx/c, 'red',label = r'$v_x /c$')
plt.plot(scale_time, vz/c, 'blue', label = r'$v_z /c$')
plt.plot(scale_time, v/c, 'orange', label = 'v/c')
plt.legend(loc = 'upper right')
matplotlib.pyplot.subplots_adjust(wspace=.75, hspace=.75)
plt.show()
plt.close()

plt.figure()
#recalculation of beta
beta_new = v / c
print('Beta : ' + str(np.mean(beta_new)))
#plot beta as a function of time
plt.plot(scale_time, beta_new)
plt.show()
plt.close()

plt.figure()
#recalculation of lorentz factor
lorentz = (1 / (2 - (2 * vz / c) - (vx ** 2) / (c ** 2)) ** .5)
#plot the lorentz factor as a function of time
plt.plot(scale_time,lorentz)
print('Lorentz factor: ' + str(np.mean(lorentz)))
#plt.savefig('C:/Users/johna/Desktop/Alpha/Velocity/Equal/' + str(alpha * 10) + '.png')
plt.show()
plt.close()