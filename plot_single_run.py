#!/bin/python

import pylab as plt
import math
import numpy as np
import shelve
from launch_single_run import method

def c_small(x):
    return x**(-0.5)

def c_large(x):
    return x**(-0.99)

x = np.arange(1, 200000, 10)

s = shelve.open("AR_single_run_shelf",'r')

# oem6 - online, exp = 0.6
# oem9 - online, exp = 0.9
# gam3 - ioem, gamma = 1 (!)

# est

plt.axis([0,200000,3,7])
plt.axhline(5.5)
plt.ylabel('estimate of sigv')
plt.plot(s['oem6'].sigv_est,'r-')
plt.plot(s['oem9'].sigv_est,'g-')
plt.plot(s['gam3'].sigv_est,'b-')
plt.savefig("sigv_est.png")
plt.close()

# eta

plt.axis([0,200000,0.00001,.1])
plt.ylabel('eta of sigv')
plt.yscale("log")
plt.plot(s['gam3'].sigv_eta,'b-')
plt.plot(x, c_small(x), 'r--')
plt.plot(x, c_large(x), 'g--')
plt.savefig("sigv_eta.png")
plt.close()

# memlen

plt.axis([0,200000,0,15000])
plt.ylabel('memlen of sigv')
plt.plot(s['gam3'].sigv_memlen,'b-')
plt.savefig("sigv_memlen.png")
plt.close()

s.close()
