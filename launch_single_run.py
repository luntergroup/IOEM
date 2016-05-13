# python script for comparing em methods over large number of runs

import random
import math
import itertools
import sys
import smc_src_single_run as smc
from operator import add, mul, sub, div
import pylab as plt
import shelve

class method:
    
    def __init__(self, em, batch_size=500, oem_exp=.75, ioem_alpha=.5, ioem_gam=1):
        self.em = em
        self.oem_exp = oem_exp
        self.batch_size = batch_size
        self.ioem_alpha = ioem_alpha
        self.ioem_gam = ioem_gam
        self.a_est = []
        self.sigw_est = []
        self.sigv_est = []
        self.dis = []
        if self.em=="ioem":
            self.a_eta = []
            self.sigw_eta = []
            self.sigv_eta = []
            self.a_memlen = []
            self.sigw_memlen = []
            self.sigv_memlen = []
	
	
    def get_info(self):
        if self.em == "batch":
            print "batch", self.batch_size
        if self.em == "online":
            print "online", self.oem_exp
        if self.em == "ioem":
            print "ioem", ioem_alpha, ioem_gam

N = 100

m1 = method("batch")
m2 = method("batch", batch_size = 1000)
m3 = method("online", oem_exp =.6)
m4 = method("online", oem_exp = .9)
m5 = method("ioem", ioem_gam=)

methods = [m1,m2,m3,m4,m5]

Start1 = smc.StartingProb("AR1", [0,5])
true_params = smc.set_params("AR1", [.95,1,5.5])
input_params = smc.set_params("AR1", [.95,3,1])


for m in methods:
    print m.em
    for run in [1]:
        random.seed(run)
        hmm = smc.HMM(200000, "AR1", Start1, true_params)
	
        record = smc.pf("AR1", obs = hmm.emission, N = N, initial_params = input_params, \
        sAR_true_params = true_params, em_method = m.em, lag=21, batch_size=m.batch_size, \
                        oem_exponent=m.oem_exp, ioem_alpha=m.ioem_alpha, ioem_gam=m.ioem_gam, sweep_indx=1)
	
        m.a_est = record.a_est
        m.sigw_est = record.sigw_est
        m.sigv_est = record.sigv_est
        if m.em=="ioem":
            m.a_eta = record.a_eta
            m.sigw_eta = record.sigw_eta
            m.sigv_eta = record.sigv_eta
            m.a_memlen = record.a_memlen
            m.sigw_memlen = record.sigw_memlen
            m.sigv_memlen = record.sigv_memlen
	    
        for i in range(len(record.a_est)):
            m.dis.append(abs(record.a_est[i]-true_params.a) + \
            abs(record.sigw_est[i]-true_params.sigw) + \
            abs(record.sigv_est[i]-true_params.sigv))

    
s = shelve.open("AR_single_run_shelf",writeback=True)
s['bem500'] = m1
s['bem1k'] = m2
s['oem6'] = m3
s['oem9'] = m4
s['gam3'] = m5
s.close()
