#FEL Equal Modulation Motion, John Anderson, Last Updated: 22/10/19

#importation of required packages used in program
import numpy as np
import matplotlib.pyplot as plt
import time
import matplotlib
import matplotlib.gridspec as grid

#formats the font for all maplotlib plots
font = {'family' : 'verdana', 'size'   : 11}
matplotlib.rc('font', **font)
start_time=time.time()  #create an initial time for when the code was run

#mass of electron (kg), charge of electron (C), magnetic field strength (T), speed of light (m/s),velocity of travelling particle (m/s)
m = 9.11e-31; q = 1.6e-19; B = .25;c = 2.99792458e8; v=0.9999998*c; alpha = 0
beta = v/c #speed of the electron as a factor of the speed of light
lorentz_factor=1/(np.sqrt(1-beta**2))
k = (np.pi * 2) / 2e-2  # wave number defined as 2pi divided by the wavelength (1/m)
au =(q * B) / (k * m * c) #undulator parameter


def x(t):   #x-position equation
    r = k*beta*c*t
    return (-au/(lorentz_factor*k*beta))*(np.cos(r)+((alpha/8)*np.sin(2*r)))

def z(t): #z-position equation
    r = k * beta * c * t
    return ((au**2)/(8*(lorentz_factor**2)*beta*k))*(np.sin(2*r)-((alpha/3)*(np.cos(3*r)-(3*np.cos(r))))-((alpha**2)/32)*np.sin(4*r))

def vx(t): #x-velocity equation
    r = k * beta * c * t
    return ((q*B)/(lorentz_factor*m*k))*(np.sin(r)-((alpha/4)*np.cos(2*r)))

def vz(t):  #z-velocity equation
    r = k * beta * c * t
    return (beta * c) + ((c*au**2)/(4*lorentz_factor**2))*(np.cos(2*r)+alpha*((np.sin(r)*np.cos(2*r))-((alpha*np.cos(4*r))/16)))

gs = grid.GridSpec(4, 1)
t1 = np.linspace(0, 1e-9, 10000) #time scale 0 to 1e-9s

#scaling factors
if alpha  == 0:
    max1 = max(t1)
    maxx = max(x(t1))
    maxz = max(z(t1))
    print(maxx)


alpha  = 0.2
scale_time = t1 / max1
scale_x = x(t1) / maxx
scale_z = z(t1) / maxz
print(z(0))

plt.figure(figsize=(6, 5))

# plot scaled x(t) against scaled time
plt.subplot(gs[0, 0], xlabel='t', ylabel='x(t)')
plt.plot(scale_time, scale_x, 'red')

# plot scaled z(t) against scaled time
plt.subplot(gs[1, 0], xlabel='t', ylabel='z(t)')
plt.plot(scale_time, scale_z, 'blue')

# plot scaled z(t) against scaled x(t)
plt.subplot(gs[2:, 0], xlabel=r'$z(t)-\beta ct$', ylabel='x(t)')
plt.plot(scale_z,scale_x, 'purple')

matplotlib.pyplot.subplots_adjust(wspace=.75, hspace=.75)
plt.show()
plt.close()

plt.figure(figsize=(5, 5))
v1 = ((vx(t1) ** 2) + (vz(t1) ** 2))**.5 #calculation of total speed

# plot scaled vx(t) against scaled time
plt.subplot(gs[0, :], xlabel='t', ylabel=r'$\beta_x$(t)')
plt.plot(scale_time, vx(t1)/c, 'red') #speed is % of speed of light

# plot scaled vz(t) against scaled time
plt.subplot(gs[1, :], xlabel='t', ylabel=r'$\beta_z$(t)')
plt.plot(scale_time, vz(t1)/c, 'blue')#speed is % of speed of light

# plot scaled v(t) against scaled time
plt.subplot(gs[2, :], xlabel='t', ylabel=r'$\beta$(t)')
plt.plot(scale_time,v1/c, 'purple')#speed is % of speed of light

#plot all velocities (% of speed of light) against scaled time
plt.subplot(gs[3, :], xlabel='t')
plt.plot(scale_time, vx(t1)/c, 'red',label = r'$v_x /c$')
plt.plot(scale_time, vz(t1)/c, 'blue', label = r'$v_z /c$')
plt.plot(t1, v1/c, 'orange', label = 'v/c')
plt.legend(loc = 'upper right')
matplotlib.pyplot.subplots_adjust(wspace=.75, hspace=.75)
plt.show()
plt.close()

plt.figure()
beta_new = v1 / c #recalculation of beta
print('Beta : ' + str(np.mean(beta_new)))
plt.plot(scale_time, beta_new)
plt.xlabel('t')
plt.ylabel(r'$\beta')
plt.show()
plt.close()

plt.figure()
#recalculation of lorentz factor
lorentz = (1 / (2 - (2 * vz(t1) / c) - (vx(t1) ** 2) / (c ** 2)) ** .5)
plt.plot(scale_time,lorentz)
print('Lorentz factor: ' + str(np.mean(lorentz)))
plt.xlabel('t')
plt.ylabel(r'$\beta')
plt.show()
plt.close()

#prints the difference of the start time to finish time, formatted to show '-time (s) -'
print("-%.2fs-" % (time.time() - start_time))


