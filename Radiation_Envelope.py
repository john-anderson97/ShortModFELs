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
n =1
w= -10
rho = 1e-3
A0 = 0+0j



#Short Modulation greater than or equal to 1.
if n >=0:
    if n ==0:
        def f1(t, y):  # x direction coupled equations
            tjz = ((au ** 2) / (8 * k * beta * beta*c * lorentz_factor ** 2)) * np.sin(2 * k * t)
            return [(np.exp(-1j * tjz * ((1 / (2 * rho)) + w))) * np.exp(1j * w * t)]
    #Equal Modulation Case
    elif n == 1:
        def f1(t, y):  # x direction coupled equations
            z1 = np.sin(2 * k * t)
            z2 = alpha/2 *(np.cos(3*k*t)-3*np.cos(k*t))
            z3 = alpha**2 /32 *np.sin(4*k*t)
            tjz = ((au ** 2) / (8 * k * beta * beta * c * lorentz_factor ** 2)) * (z1-z2-z3)
            return [(np.exp(-1j * tjz * ((1 / (2 * rho)) + w))) * np.exp(1j * w * t)]
    #Half Modulation Case
    elif n==2:
        def f1(t, y):  # x direction coupled equations
            z1 = np.sin(2 * k * t)
            z2 = alpha / 12 * (np.cos(4 * k * t) - 2 * np.cos(k * t))
            z3 = alpha/ 2  * np.sin(2 * k * t)
            z4 = alpha**2 / 216 *np.sin(6*k*t)
            z5 = alpha**2 / 12 *np.sin(2*k*t)
            z6 = alpha**2 / 24 *np.sin(4*k*t)
            z7 = alpha**2 / 8 *np.sin(2*k*t)
            tjz = ((au ** 2) / (8 * k * beta * beta * c * lorentz_factor ** 2)) * (z1 - z2 - z3-z4-z5-z6-z7)
            return [(np.exp(-1j * tjz * ((1 / (2 * rho)) + w))) * np.exp(1j * w * t)]
    #General Modulation Case
    else:
        def f1(t, y):  # x direction coupled equations
            z1 = np.sin(2 * k * t)
            z2 = (-2 * alpha / (n + 1)) * ((np.cos((n + 2) * k * t) / (n + 2)) + (np.cos(n * k * t) / n))
            z3 = (-2 * alpha / (n - 1)) * ((np.cos(n * k * t) / n) - (np.cos((n - 2) * k * t) / (n - 2)))
            z4 = -(alpha ** 2 / 4) * ((np.sin(2 * (n + 1) * k * t)) / ((n + 1) ** 3))
            z5 = (-alpha ** 2 / (2 * (n + 1) * (n - 1))) * (np.sin(2 * k * t) + (np.sin(2 * n * k * t) / n))
            z6 = -(alpha ** 2 / 4) * ((np.sin(2 * (n - 1) * k * t)) / ((n - 1) ** 3))
            tjz = ((au ** 2) / (8 * k * beta * beta*c * lorentz_factor ** 2)) * (z1 + z2 + z3 + z4 + z5 + z6)
            return [(np.exp(-1j * tjz * ((1 / (2 * rho)) + w))) * np.exp(1j * w * t)]
else:
    print('n must have a value greater than or equal to 1')
    sys.exit()

#RK45 method plus plotting
def main():
    tTime = 8*np.pi # total zbar
    stepsize = tTime / 1000 # number of steps

    ##solutions to equations using RK45 method with initial condidtion, timeframe and step size.
    x_sol = solve_ivp(f1, t_span=[0, tTime], y0=[A0], method='RK45', max_step=stepsize, atol=1, rtol=1)

    t1 = x_sol.t  # data for zbar
    x = x_sol.y[0]  # data for |A|

    #append to arrays
    z.append(t1)
    a.append(abs(x/(np.exp(1j*w*t1))))


a= []
z =[]
w1 = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
print(len(z))

counter = 0
while counter <90:
    main()
    print(counter)
    i = 0
    while i < len(a[counter]):
        w1[counter].append(w)
        i=i+1
    i=0
    w = w+(20/90)
    counter = counter +1



#plot of omegabar and zbar against |A|
fig = plt.figure(1)
sc = plt.scatter(z,w1, 10, a)
plt.xlabel(r'$\bar{z}$')
plt.ylabel(r'$\bar{\omega}$')
cbar = plt.colorbar(sc)
cbar.set_label(r'$|A|$')
plt.show()

