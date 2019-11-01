# python script for comparing em methods over large number of runs

import random
import math
import itertools
import sys
import smc_src as smc
from operator import add, mul, sub, div
import pylab as plt
import shelve

class method:
    
    def __init__(self, em, batch_size=500, oem_exp=.75, ioem_gam=1, ioem_alpha = None, min_mem = 500):
        self.em = em
        self.oem_exp = oem_exp
        self.batch_size = batch_size
        self.ioem_alpha = ioem_alpha
        self.ioem_gam = ioem_gam
        self.a_est = []
        self.sigw_est = []
        self.sigv_est = []
        self.dis = []
        self.min_mem = min_mem
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


if __name__ == "__main__":

    N = 5000

    ## read observations from a column of the datafile
    infile = open( sys.argv[1], "r")    ## argument: gbpusd-logreturns.txt
    infile.readline()
    incol = 1 + int(sys.argv[2])        ## argument: 0-23, refers to hour column
    obs = [ float(line.split('\t')[ incol ]) for line in infile ]

    gammas = [ 0.5, 1, 1.5, 2, 3 ]
    methods = [ method("ioem", ioem_gam = g, min_mem = 100) for g in gammas ]    ## scaling of lookback.
    methodnames = [ 'ioem{}'.format(g) for g in gammas ]

    Start1 = smc.StartingProb("AR1", [0,0.005])      ## not used
    true_params = smc.set_params("AR1", [.95,1,5.5]) ## not used
    input_params = smc.set_params("AR1", [0.05, 0.005, 0.005])

    for m in methods:
        print m.em
        for run in [1]:
            random.seed(run)

            #hmm = smc.HMM(200000, "AR1", Start1, true_params)
            #obs = hmm.emission

            record = smc.pf("AR1", obs = obs, N = N, initial_params = input_params, \
                                sAR_true_params = true_params, em_method = m.em, lag=21, batch_size=m.batch_size, \
                                oem_exponent=m.oem_exp, ioem_gam=m.ioem_gam, sweep_indx=1, min_mem = m.min_mem)

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

    
    s = shelve.open("AR_single_run_shelf_col{}".format(sys.argv[2]),writeback=True)
    for name, method in zip(methodnames, methods):
        s[ name ] = method
    s.close()
