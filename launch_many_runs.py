# python script for comparing em methods over large number of runs

import random
import math
import itertools
import sys
import smc_src_many_runs as smc
from operator import add, mul, sub, div
import pylab as plt
import shelve

class method:
    
    def __init__(self, em, batch_size=500, oem_exp=.75, ioem_gam=1):
        self.em = em
        self.oem_exp = oem_exp
        self.batch_size = batch_size
        self.ioem_alpha = ioem_alpha
        self.ioem_gam = ioem_gam
        self.a1k_estimates = []
        self.sigw1k_estimates = []
        self.sigv1k_estimates = []
        self.dis1k = []
        self.a5k_estimates = []
        self.sigw5k_estimates = []
        self.sigv5k_estimates = []
        self.dis5k = []
        self.a10k_estimates = []
        self.sigw10k_estimates = []
        self.sigv10k_estimates = []
        self.dis10k = []
        self.a20k_estimates = []
        self.sigw20k_estimates = []
        self.sigv20k_estimates = []
        self.dis20k = []
        self.a50k_estimates = []
        self.sigw50k_estimates = []
        self.sigv50k_estimates = []
        self.dis50k = []
	
    def get_info(self):
        if self.em == "batch":
            print "batch", self.batch_size
        if self.em == "online":
            print "online", self.oem_exp
        if self.em == "ioem":
            print "ioem", ioem_alpha, ioem_gam

N = 100


m1 = method("batch")
m2 = method("batch", batch_size = 5000)
m3 = method("online", oem_exp =.6)
m4 = method("online", oem_exp = .9)
m5 = method("ioem")

methods = [m1,m2,m3,m4,m5]

Start1 = smc.StartingProb("AR1", [0,5])
true_params = smc.set_params("AR1", [.95,1,5.5])
input_params = smc.set_params("AR1", [.8,3,1])


for m in methods:
    print m.em
    for run in range(100):
        random.seed(run)
        if run%10==0:
            print run
        hmm = smc.HMM(50000, "AR1", Start1, true_params)
        
        end_pars = smc.pf("AR1", obs = hmm.emission, N = N, initial_params = input_params, \
            sAR_true_params = None, em_method = m.em, lag=21, batch_size=m.batch_size, \
            oem_exponent=m.oem_exp, ioem_gam=m.ioem_gam, sweep_indx=1)
        
        m.a1k_estimates.append(end_pars[0].a.estimate)
        m.sigw1k_estimates.append(end_pars[0].sigw.estimate)
        m.sigv1k_estimates.append(end_pars[0].sigv.estimate)
        m.dis1k.append(abs(end_pars[0].a.estimate-true_params.a)+ \
        abs(end_pars[0].sigw.estimate-true_params.sigw)+ \
        abs(end_pars[0].sigv.estimate-true_params.sigv))
        
        m.a5k_estimates.append(end_pars[1].a.estimate)
        m.sigw5k_estimates.append(end_pars[1].sigw.estimate)
        m.sigv5k_estimates.append(end_pars[1].sigv.estimate)
        m.dis5k.append(abs(end_pars[1].a.estimate-true_params.a)+ \
        abs(end_pars[1].sigw.estimate-true_params.sigw)+ \
        abs(end_pars[1].sigv.estimate-true_params.sigv))
        
        m.a10k_estimates.append(end_pars[2].a.estimate)
        m.sigw10k_estimates.append(end_pars[2].sigw.estimate)
        m.sigv10k_estimates.append(end_pars[2].sigv.estimate)
        m.dis10k.append(abs(end_pars[2].a.estimate-true_params.a)+ \
        abs(end_pars[2].sigw.estimate-true_params.sigw)+ \
        abs(end_pars[2].sigv.estimate-true_params.sigv))
        
        m.a20k_estimates.append(end_pars[3].a.estimate)
        m.sigw20k_estimates.append(end_pars[3].sigw.estimate)
        m.sigv20k_estimates.append(end_pars[3].sigv.estimate)
        m.dis20k.append(abs(end_pars[3].a.estimate-true_params.a)+ \
        abs(end_pars[3].sigw.estimate-true_params.sigw)+ \
        abs(end_pars[3].sigv.estimate-true_params.sigv))
        
        m.a50k_estimates.append(end_pars[4].a.estimate)
        m.sigw50k_estimates.append(end_pars[4].sigw.estimate)
        m.sigv50k_estimates.append(end_pars[4].sigv.estimate)
        m.dis50k.append(abs(end_pars[4].a.estimate-true_params.a)+ \
        abs(end_pars[4].sigw.estimate-true_params.sigw)+ \
        abs(end_pars[4].sigv.estimate-true_params.sigv))	

    
s = shelve.open("AR_many_runs_shelf",writeback=True)
s['bem500'] = m1
s['bem5k'] = m2
s['oem6'] = m3
s['oem9'] = m4
s['ioem'] = m5
s.close()
