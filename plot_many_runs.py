#!/bin/python

import pylab as plt
import math
import shelve
from launch_many_runs import method

#sigv

plt.rcParams.update({'font.size': 20})
plt.rcParams['xtick.major.pad']='12'
plt.axis([0,5,3,7])
plt.axhline(5.5)
data = []

for m in ['bem500','bem5000','oem0.6','oem0.9','ioem']:
    s = shelve.open("AR_many_runs_shelf",flag='r')
    data.append(s[m].__dict__['sigv50k_estimates'])
    s.close()
plt.boxplot(data)

plt.axvline(2.5,color="grey")
plt.axvline(4.5,color="grey")
plt.gcf().subplots_adjust(left=0.12, bottom=0.2)
plt.ylabel(r"$\hat{\sigma}_{v,50k}$",rotation=90,fontsize=24)

plt.xticks([1,2,3,4], [r"$500$",r"$5k$",r"$0.6$",r"$0.9$"],fontsize=20,rotation=90)
plt.xlabel("tuning parameter value", fontsize=20)

plt.tick_params(axis='x', which='both', length=0)

plt.twiny()
plt.xlim(0,5)
plt.xticks([1,3,4.5], ["BEM","OEM","IOEM"], fontsize=20)
plt.tick_params(axis='x', which='both', length=0)

plt.savefig("sigv_50k.png")
plt.close()

#sigw

plt.rcParams.update({'font.size': 20})
plt.rcParams['xtick.major.pad']='12'
plt.axis([0,5,0,5])
plt.axhline(1)
data = []

for m in ['bem500''bem5000','oem0.6','oem0.9','ioem']:
    s = shelve.open("AR_many_runs_shelf",flag='r')
    data.append(s[m].__dict__['sigw50k_estimates'])
    s.close()
plt.boxplot(data)

plt.axvline(2.5,color="grey")
plt.axvline(4.5,color="grey")
plt.gcf().subplots_adjust(left=0.12, bottom=0.2)
plt.ylabel(r"$\hat{\sigma}_{w,50k}$",rotation=90,fontsize=24)

plt.xticks([1,2,3,4], [r"$500$",r"$5k$",r"$0.6$",r"$0.9$"],fontsize=20,rotation=90)
plt.xlabel("tuning parameter value", fontsize=20)

plt.tick_params(axis='x', which='both', length=0)

plt.twiny()
plt.xlim(0,5)
plt.xticks([1,3,4.5], ["BEM","OEM","IOEM"], fontsize=20)
plt.tick_params(axis='x', which='both', length=0)

plt.savefig("sigw_50k.png")
plt.close()

# a

plt.rcParams.update({'font.size': 20})
plt.rcParams['xtick.major.pad']='12'
plt.axis([0,5,.3,1])
plt.axhline(.95)
data = []

for m in ['bem500','bem5000','oem0.6','oem0.9','ioem']:
    s = shelve.open("AR_many_runs_shelf",flag='r')
    data.append(s[m].__dict__['a50k_estimates'])
    s.close()
plt.boxplot(data)

plt.axvline(2.5,color="grey")
plt.axvline(4.5,color="grey")
plt.gcf().subplots_adjust(left=0.12, bottom=0.2)
plt.ylabel(r"$\hat{a}_{50k}$",rotation=90,fontsize=24)

plt.xticks([1,2,3,4], [r"$500$",r"$5k$",r"$0.6$",r"$0.9$"],fontsize=20,rotation=90)
plt.xlabel("tuning parameter value", fontsize=20)

plt.tick_params(axis='x', which='both', length=0)

plt.twiny()
plt.xlim(0,5)
plt.xticks([1,3,4.5], ["BEM","OEM","IOEM"], fontsize=20)
plt.tick_params(axis='x', which='both', length=0)

plt.savefig("a_50k.png")
plt.close()
