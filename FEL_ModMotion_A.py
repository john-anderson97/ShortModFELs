#FEL Modulation Motion, John Anderson, Last Updated: 22/10/19

#importation of required packages used in program
import numpy as np
import matplotlib.pyplot as plt
import time
import matplotlib
import matplotlib.gridspec as gridspec

#formats the font for all maplotlib plots
font = {'family' : 'verdana', 'size'   : 11}
matplotlib.rc('font', **font)

start_time=time.time()  #create an initial time for when the code was run

#mass of electron (kg), charge of electron (C), magnetic field strength (T), speed of light (m/s),velocity of travelling particle (m/s), amplitude of variation
m = 9.11e-31; q = 1.6e-19; B = .25;c = 2.99792458e8; v=0.9999998*c; alpha = 0
beta = v/c #speed of the electron as a factor of the speed of light
lorentz_factor=1/(np.sqrt(1-beta**2))
k = (np.pi * 2) / 2e-2  # wave number defined as 2pi divided by the wavelength (1/m)
au =(q * B) / (k * m * c)
km = k*.25

def x(t):
    x1 = np.sin((k+km)*beta*c*t)/((k+km)**2)
    x2 = np.sin((k-km)*beta*c*t)/((k-km)**2)
    x3 = (alpha*k*au)/(2*lorentz_factor*beta)
    return (-au/(lorentz_factor*k*beta))*np.cos(k*beta*c*t)-x3*(x1+x2)

def z(t):
    x1 = np.sin(2*(km+k)*beta*c*t)/((k+km)**3)
    x2 = np.sin(2 * (km - k) * beta * c * t) / ((k - km) ** 3)
    x3 = (1/(32*beta))*((alpha*k*au)/(lorentz_factor))**2
    return -(((au**2)/(8*(lorentz_factor**2)*beta*k))*np.sin(2*k*beta*c*t))+x3*(x1+x2)
##### speed equations need updated
def vx(t):
    x1 = (q*B)/(lorentz_factor*m)
    x2 = np.sin(k*beta*c*t)/k
    x3 = alpha/2
    x4 = np.cos((km+k)*beta*c*t)/(km+k)
    x5 = np.cos((km-k)*beta*c*t)/(km-k)
    return x1*(x2-(x3*(x4+x5)))

def vz(t):    #z direction coupled equations
    x1 = (c*au**2)/(4*lorentz_factor**2)
    x2 = np.cos(2*k*beta*c*t)
    x3 = (c*(alpha**2)*(k**2)*(au**2))/(16*lorentz_factor**2)
    x4 = np.cos(2*(km+k)*beta*c*t)/((km+k)**2)
    x5 = np.cos(2*(km-k)*beta*c*t)/((km-k)**2)
    return beta*c +(x1*x2)-(x3*(x4-x5))

#arrays for comparision of max amplitudes of each alpha
x_array = []
z_array= []
lor_array =[]
alpha_array =[]
beta_array = []

#loop until alpha breaks the system, alpha must be less than 1
while alpha <=1:
    gs = gridspec.GridSpec(4, 1)
    t1 = np.linspace(0, 1e-8, 10000)

    if alpha == 0:
        #scaling values using maximum value when alpha = 0
        scalex = max(x(t1))
        scaletime = max(t1)
        scalez = max(z(t1))

    x_array.append(max(x(t1) / scalex))
    z_array.append(max(z(t1) / scalez))
    alpha_array.append(alpha)

    plt.figure(figsize=(6, 5))

    # plot scaled x(t) against scaled time
    plt.subplot(gs[0, 0], xlabel='t', ylabel='x(t)')
    plt.plot(t1 / scaletime, x(t1) / scalex, 'red')

    # plot scaled z(t) against scaled time
    plt.subplot(gs[1, 0], xlabel='t', ylabel='z(t)')
    plt.plot(t1 / scaletime, z(t1) / scalez, 'blue')

    # plot scaled z(t) against scaled x(t)
    plt.subplot(gs[2:, 0], xlabel=r'$z(t)-\beta ct$', ylabel='x(t)')
    plt.plot(z(t1) / scalez, x(t1) / scalex, 'purple')

    matplotlib.pyplot.subplots_adjust(wspace=.75, hspace=.75)

    plt.suptitle(r'$\alpha = $' + str(alpha))
    plt.savefig('C:/Users/johna/Desktop/Alpha/Position/Equal/' + str(alpha * 10) + '.png')
    #plt.show()
    plt.close()

    plt.figure(figsize=(10, 5))
    #plot of x-velocity against time

    #plt.plot(t1/scaletime, vx(t1)/c, 'red',label = r'$v_x /c$')
    # plot of y-velocity against time
    #plt.plot(t1/scaletime, vz(t1)/c, 'blue', label = r'$v_z /c$')

    #calculation of total speed
    v = (vx(t1) ** 2) + (vz(t1) ** 2)

    # plot of velocity against time
    plt.plot(t1, vz(t1), 'orange', label = 'v/c')
    plt.legend(loc = 'upper right')
    #recalculation of beta
    #plt.show()
    beta_new = v / c
    beta_array.append(np.mean(beta_new))

    #recalculation of lorentz factor


    beta_2 = (v ** .5) / c
    print(beta_2)
    lorentz = ((1 / (1 - (beta_2 ** 2)) ** .5))
    lor_array.append(np.mean(lorentz))
    print(lorentz)

    #plt.savefig('C:/Users/johna/Desktop/Alpha/Velocity/Equal/' + str(alpha * 10) + '.png')
    #plt.show()
    plt.close()

    #increase alpha by .1
    alpha = round(alpha + .1, 2)

plt.figure()
#plot of max x and max z for each alpha, scaled according to the alpha = 0
plt.plot(alpha_array, x_array, label = 'Max x')
plt.plot(alpha_array, z_array, label = 'Max z')
plt.legend(loc = 'upper left')
plt.xlabel(r'$\alpha$')
plt.show()

plt.figure()
#plot of lorentz factor for each alpha
plt.plot(alpha_array, lor_array)
plt.xlabel(r'$\alpha$')
plt.ylabel(r'$\gamma$')
plt.show()

plt.figure()
#plot of beta for each alpha
plt.plot(alpha_array, beta_array)
plt.ylim(0.99,1.01)
plt.xlabel(r'$\alpha$')
plt.ylabel(r'$\beta$')
plt.show()
#prints the difference of the start time to finish time, formatted to show '-time (s) -'
print("-%.2fs-" % (time.time() - start_time))

#shows the plots in the window
## Time to run code: 0.28s