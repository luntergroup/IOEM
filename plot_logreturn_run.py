#!/bin/python

import pylab as plt
import math
import numpy as np
import shelve
from launch_logreturn_run import method

##
## read observations (00:00)
##
infile = open( "gbpusd-logreturns.txt", "r" )
infile.readline()
incol = 1 + 0
obs, gbpusd = zip( *[ (float(line.split('\t')[ incol ]),
                       float(line.split('\t')[ 25 ]) ) for line in infile ] )


##
## from the ARMA(1,1) inference run in an R session
##
arma_relcov = [-0.00556616511471543,-0.0100981576683794,-0.0131099141220133,-0.0219624786546186,-0.0251946437224198,-0.0217877122208273,-0.0167980352406851,-0.0102731996344379,-0.00856463497124997,-0.0077853136486904,0.00450804989078992,0.00278569953584888,-0.0165393606404879,-0.0237926656380352,-0.00997867981624065,-0.00405634424446305,-0.012847380370993,-0.0338358144293583,-0.039571774951013,-0.0201024810808082,-0.0168759608217075,-0.018985090991212,-0.0219938096255779,-0.000921611501839996]


##
## read IOEM inferences
##
s = [ shelve.open("AR_single_run_shelf_col{}".format(c))
      for c in range(24) ]
gammas = [0.5, 1, 1.5, 2, 3]
methodnames = [ 'ioem{}'.format(g) for g in gammas ]
methodidx = 3

##
## nice colour sequence
##
colormap = plt.get_cmap('jet')
custom_cycler = (plt.cycler(color = [ colormap(c/23.0) for c in range(0,24) ]) +
                 plt.cycler(linestyle = ['-','-.',':'] * 8 ))
plt.rc('axes', prop_cycle = custom_cycler)


##
## plot data, and online inferences of the 3 parameters
##
plt.figure(figsize=(10,3))
ax = plt.subplot(1, 4, 1)
ax.locator_params(nbins=5, axis='y')
plt.yticks(rotation=90)
plt.axis([0, 1499, -0.02, 0.02])
#plt.ylabel("GBP/USD log day return")
plt.title("GBP/USD log day returns", fontsize=10)
plt.xlabel("time (days)")
plt.plot( obs, 'o' , markeredgecolor="grey", markersize = 0.5, linestyle="-", linewidth=0.3, color="0.75" )

ax2 = ax.twinx()
ax2.locator_params(nbins=5, axis='y')
ax2.tick_params(axis='y', rotation=90)
ax2.plot( gbpusd, linestyle="-", linewidth=1, color=[0.3,0.3,0.7,0.5] )


yy = np.transpose( [ s[i][ methodnames[methodidx] ].a_est for i in range(24) ] )
ax = plt.subplot(1, 4, 2)
ax.locator_params(nbins=4, axis='y')
plt.yticks(rotation=90)
plt.axis([0, 1499, -0.26, 0.1])
plt.xlabel("time (days)")
plt.title("$\hat{a}_t$", fontsize=10)
plt.plot(yy, linewidth=0.3)
plt.axhline(0, color="gray", dashes=[0.5,0.5], linewidth=1)


yy = np.transpose( [ s[i][ methodnames[methodidx] ].sigw_est for i in range(24) ] )
ax = plt.subplot(1, 4, 3)
ax.locator_params(nbins=4, axis='y')
plt.yticks(rotation=90)
plt.axis([0, 1499, 0.001, 0.006])
plt.title("$\hat{\sigma}_{w,t}$", fontsize=10)
plt.xlabel("time (days)")
plt.plot(yy, linewidth=0.3)


yy = np.transpose( [ s[i][ methodnames[methodidx] ].sigv_est for i in range(24) ] )
ax = plt.subplot(1, 4, 4)
ax.locator_params(nbins=4, axis='y')
plt.yticks(rotation=90)
plt.axis([0, 1499, 0.001, 0.006])
plt.xlabel("time (days)")
plt.title("$\hat{\sigma}_{v,t}$", fontsize=10)
plt.plot(yy, linewidth=0.3)

plt.subplots_adjust(top=0.92, bottom=0.15, left=0.04, right=0.98, hspace=0.0,
                    wspace=0.25)

plt.savefig("usdgbp.png")
plt.close()

#
# calculate cov/var of the 24 daily return sequences
#

rel_cov = []
for i in range(24):
    model = s[i][ methodnames[methodidx] ]
    a, sigw, sigv = model.a_est[-1], model.sigw_est[-1], model.sigv_est[-1]
    var = sigv + sigw / (1 - a*a)
    cov = sigw * a / (1 - a*a)
    rel_cov.append( cov/var )
#print "c({})".format( ",".join(map(str, rel_cov)))


plt.axis([-0.042, 0.006, min(rel_cov)-0.005, max(rel_cov)+0.005])
plt.ylabel("correlation - IOEM")
plt.xlabel("correlation - ARMA(1,1)")
plt.plot( arma_relcov, rel_cov, 'b.' )
plt.axhline( 0, color="gray", dashes=[0.5,0.5], linewidth=2)
plt.axvline( 0, color="gray", dashes=[0.5,0.5], linewidth=2)
plt.savefig("correlation.png")
plt.close()

#
# memlen
#

plt.figure(figsize=(10,3))
vars = ["a_memlen", "sigw_memlen", "sigv_memlen"]

for idx, var in enumerate(vars):
    yy = eval( "np.transpose( [ s[i][ methodnames[methodidx] ].{} for i in range(24) ] )".format(var) )
    ax = plt.subplot(1, 3, idx+1)
    ax.locator_params(nbins=4, axis='y')
    plt.yticks(rotation=90)
    plt.axis([0,1500,0,750])
    plt.ylabel('{}'.format(var) )
    plt.plot( yy, linewidth=0.3 )

plt.savefig("memlen.png")
plt.close()




