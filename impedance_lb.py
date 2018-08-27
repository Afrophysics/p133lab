'''This code is dedicated to the DC circuit observations for a battery with an uknown internal resistance. The goal of this code 
is to calculate the potential within the battery and the resistance of the internal resistor using a linear fit model for the IV
diagram of this'''

import numpy as np
import matplotlib.pyplot as plt

#Set up recorded parameters for applied external resistance from variable resistor box and 
#the observed potential for the given resistance
R = np.array([10.2, 40.1, 200, 990, 40000])
err_r = np.repeat(0.3, len(R))
V = np.array([0.111, 0.420, 1.760, 4.91, 8.72])
err_v = np.array([0.001, 0.001, 0.001, 0.01, 0.01])
I = V/R
err_cur = ((err_v**2 + (I * err_r)**2)**0.5)/R




print("Part (a) is meant to build the algorithm")

#Expectation value of the nth-power of the x-array weighted by the inverse of the y-variance
def u(x, sig, n):
    tp = x**n
    bt = sig**2
    ary = tp/bt
    return np.sum(ary)

#Expectation value of the inner product between the the nth-power of the x-array and the y-array weighted by the inverse of the y-variance
def w(x, y, sig, n):
    tp = y*(x**n)
    bt = sig**2
    ary = tp/bt
    return np.sum(ary)
    
#Linear-fit algorithm using least-squares method
def linfit(x_data, y_data, y_error):
    u_0 = u(x_data, y_error, 0)
    u_1 = u(x_data, y_error, 1)
    u_2 = u(x_data, y_error, 2)
    
    w_0 = w(x_data, y_data, y_error, 0)
    w_1 = w(x_data, y_data, y_error, 1)
    w_2 = w(x_data, y_data, y_error, 2)
    
    D = (u_0 * u_2) - (u_1 ** 2)
    stnd_1 = (u_0 * w_1) - (w_0 * u_1)
    stnd_2 = (u_2 * w_0) - (w_1 * u_1)
    
    slope = stnd_1/D
    var_sl = u_0/D
    
    y_intr = stnd_2/D
    var_yt = u_2/D
    
    return np.array([slope, y_intr, var_sl, var_yt])
    
print("Part (a) complete")
    
    
print("Part (b) is shall test the algorithm using the sample data below:")

x = I
y = V
yerr = err_v

print('x data points: %s'%x)
print('y data points: %s'%y)
print('y-error data points: %s'%yerr)

rs = linfit(x, y, yerr)
print("results")
print("slope = %.3f +/- %.3f"%(rs[0], rs[2]**0.5))
print("intercept = %.3f +/- %.3f"%(rs[1], rs[3]**0.5))

print("Part (b) complete")
print('')
print("Part (c) we shall plot our results")
import matplotlib.pyplot as plt

fig = plt.figure()
plt.errorbar(x, y, yerr = yerr, fmt='.', label='measured volts')
smpl_x = np.linspace(0, 0.0109, 31)

#plot the linear approximation
f_x = rs[1] + (rs[0] * smpl_x)
plt.plot(smpl_x, f_x, '--', label='%.2f + %.2f I'%(rs[1], rs[0]))
plt.grid(True)
plt.title('DC Lab')
plt.legend()
plt.xlabel('Current')
plt.ylabel('Potential')
plt.show()
#fig.savefig('Linear_fit_ex.png')

print("Problem 2 complete")
