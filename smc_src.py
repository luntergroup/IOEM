# want this to be general enough to easily translate to new models

import random
import math
import copy
from operator import add, mul, sub, div


### hmm functions

def Start(model, StartingProbs):
    '''
    create X_1 for the given model with the specified initial distribution (mu)
    '''

    if model=="coin":
        r = random.random()
        if(StartingProbs.f > r):
            x1 = 'F'
        else:
            x1 = 'U'

    if model=="AR1":
        x1 = random.gauss(StartingProbs.mu, StartingProbs.sigma)

    if model=="sAR1" or model=="sar":
        x1 = random.gauss(StartingProbs.mu, StartingProbs.sigma)

    if model=="sv":
        # want N(0, sqrt(sigma2/(1-phi2))), set this in StartingProbs
        x1 = random.gauss(StartingProbs.mu, math.sqrt(StartingProbs.var))

    if model=="two_ar":
        x1 = [random.gauss(StartingProbs.mu1, StartingProbs.sigma1),
        random.gauss(StartingProbs.mu2, StartingProbs.sigma2)]

    if model=="two_ar_shared_sigv":
        x1 = [random.gauss(StartingProbs.mu1, StartingProbs.sigma),
        random.gauss(StartingProbs.mu2, StartingProbs.sigma)]

    if model=="ar_sigw":
        x1 = random.gauss(StartingProbs.mu, StartingProbs.sigma)

    return x1

def Transition(model, params, old_x, sAR_true_params):
    '''
    create X_{t+1} from X_{t}
    sAR_true_params is only used in the sAR model, but necessary to run any as currently written
    '''
# sAR uses the true parameters for a, sigw as those are known in that model.

    if model=="coin":
        r = random.random()
        if old_x=='F':
            if(params.ff > r):
                new_x='F'
            else:
                new_x='U'
        else:
            if(params.uf > r):
                new_x='F'
            else:
                new_x='U'

    if model=="AR1":
        new_x = random.gauss(params.a * old_x, params.sigw)

    if model=="sAR1" or model=="sar":
        new_x = random.gauss(sAR_true_params.a * old_x, sAR_true_params.sigw)

    if model=="sv":
        new_x = random.gauss(params.phi*old_x, params.sig)

    if model=="two_ar" or model=="two_ar_shared_sigv":
        new_x = [random.gauss(params.a1 * old_x[0], params.sigw1),
        random.gauss(params.a2 * old_x[1], params.sigw2)]

    if model=="ar_sigw":
        new_x = random.gauss(old_x, params.sigw)

    return new_x

def Emission(model, params, x):
# create Y_t

    if model=="coin":
        r = random.random()
        if x=='F':
            if(params.ft > r):
                y='T'
            else:
                y='H'
        else:
            if(params.ut > r):
                y='T'
            else:
                y='H'

    if model=="AR1":
        y=random.gauss(x, params.sigv)

    if model=="sAR1":
        y=random.gauss(x, params.sigv)

    if model=="sar":
        y=random.gauss(x, params.sigv2**(.5))

    if model=="sv":
        y=random.gauss(0, params.b * math.exp(x/2.0))

    if model=="two_ar":
        y=[random.gauss(x[0], params.sigv1),
        random.gauss(x[1],params.sigv2)]

    if model=="two_ar_shared_sigv":
        y=[random.gauss(x[0], params.sigv),
        random.gauss(x[1],params.sigv)]

    if model=="ar_sigw":
        y=random.gauss(x, 0.5)

    return y

class HMM:
# contains true X,Y data

    def __init__(self, nSim, model, StartingProb, true_params):
        self.model = model
        self.states = []
        self.emission = []
        self.states.append(Start(model, StartingProb))
        for i in range(1,nSim):
            self.states.append(Transition(model, true_params, self.states[-1], true_params))
        for i in range(nSim):
            self.emission.append(Emission(model, true_params, self.states[i]))

    def get_model(self):
        print self.model


class particle:
# contains state, unnormalized weight, relevant state history of a particle

    def __init__(self, model, N, params, obs ):

        self.weight = 1.0/N

        if model == "coin":
            # sample x|y using initial emission params
            r = random.random()

            if obs == 'T':              #TAIL
                pFgivenT = params.ft.estimate/(params.ft.estimate+params.ut.estimate)
                if pFgivenT > r:
                    self.state = 'F'
                else:
                    self.state = 'U'

            else:                       #HEAD
                pFgivenH = (1-params.ft.estimate)/(2-params.ft.estimate-params.ut.estimate)
                if pFgivenH > r:
                    self.state = 'F'
                else:
                    self.state = 'U'

        if model == "AR1":
            # sample x|y using initial emission params
            self.state = random.gauss(obs,params.sigv.estimate)

        if model == "sAR1":
            # sample x|y using initial emission params
            self.state = random.gauss(obs,params.sigv.estimate)

        if model == "sar":
            self.state = random.gauss(obs,params.sigv2.estimate**(.5))

        if model == "sv":
            #initialize without data:
            self.state = random.gauss(0,params.sig.estimate**2/(1-params.phi.estimate**2))
            #initialize non 1/N weights:
            self.weight = ( ((2*math.pi*params.b.estimate**2*math.exp(self.state))**(-.5)) * \
            math.exp(-.5 * (obs**2/(params.b.estimate*math.exp(self.state/2.0))**2)) )
            #check these are normalized!!!!!

        if model == "two_ar":
            # sample x|y using initial emission params
            self.state = [random.gauss(obs[0],params.sigv1.estimate),
            random.gauss(obs[1],params.sigv2.estimate)]

        if model == "two_ar_shared_sigv":
            # sample x|y using initial emission params
            self.state = [random.gauss(obs[0],params.sigv.estimate),
            random.gauss(obs[1],params.sigv.estimate)]

        if model == "ar_sigw":
            self.state = random.gauss(obs,0.5)


        self.state_hist = [self.state]

    def update_particle(self, model, params, obs, lag, sAR_true_params):
        # sample x_t|x_t-1 and update particle weight

        if model == "coin":
            r = random.random()

            if self.state == 'F':      #Fair coin ->
                if params.ff.estimate > r: # -> Fair coin
                    self.state = 'F'
                    # update weight dependent upon observation
                    if obs == 'T':     #TAIL
                        self.weight = self.weight * params.ft.estimate
                    else:            #HEAD
                        self.weight = self.weight * (1-params.ft.estimate)
                else:
                    self.state = 'U'
                    # update weight dependent upon observation
                    if obs == 'T':     #TAIL
                        self.weight = self.weight * params.ut.estimate
                    else:            #HEAD
                        self.weight = self.weight * (1-params.ut.estimate)

            else:                    #Unfair coin ->
                if params.uf.estimate > r:
                    self.state = 'F'
                    # update weight dependent upon observation
                    if obs == 'T':     #TAIL
                        self.weight = self.weight * params.ft.estimate
                    else:            #HEAD
                        self.weight = self.weight * (1-params.ft.estimate)
                else:
                    self.state = 'U'
                    # update weight dependent upon observation
                    if obs == 'T':     #TAIL
                        self.weight = self.weight * params.ut.estimate
                    else:            #HEAD
                        self.weight = self.weight * (1-params.ut.estimate)

        if model == "AR1":
            self.state = random.gauss(params.a.estimate * self.state, params.sigw.estimate)
            self.weight = self.weight * (2*math.pi*(params.sigv.estimate**2))**(-.5) * \
             math.exp(-.5 * (self.state - obs)**2/params.sigv.estimate**2)

        if model == "sAR1":
            self.state = random.gauss(sAR_true_params.a * self.state, sAR_true_params.sigw)
            self.weight = self.weight * (2*math.pi*(params.sigv.estimate**2))**(-.5) * \
             math.exp(-.5 * (self.state - obs)**2/params.sigv.estimate**2)

        if model == "sar":
            self.state = random.gauss(sAR_true_params.a * self.state, sAR_true_params.sigw)
            self.weight = self.weight * (2*math.pi*(params.sigv2.estimate))**(-.5) * \
             math.exp(-.5 * (self.state - obs)**2/params.sigv2.estimate)

        if model == "sv":
            self.state = random.gauss(params.phi.estimate*self.state, params.sig.estimate)
            #checked:
            self.weight = self.weight * ( ((2*math.pi*params.b.estimate**2*math.exp(self.state))**(-.5)) \
            * math.exp(-.5 * (obs**2/(params.b.estimate*math.exp(self.state/2.0))**2)) )

        if model == "two_ar":
            self.state = [random.gauss(params.a1.estimate * self.state[0], params.sigw1.estimate), \
                random.gauss(params.a2.estimate * self.state[1], params.sigw2.estimate)]
            # below is proportional to the density of two independent gaussians
            # always normalize weights relative to each other, so no need to normalize here
            # removed 1/sqrt(2pi)
            self.weight = self.weight * (params.sigv1.estimate**(-1)) * \
             math.exp(-.5 * (self.state[0] - obs[0])**2/params.sigv1.estimate**2) * \
             (params.sigv2.estimate**(-1)) * \
             math.exp(-.5 * (self.state[1] - obs[1])**2/params.sigv2.estimate**2)

        if model == "two_ar_shared_sigv":
            self.state = [random.gauss(params.a1.estimate * self.state[0], params.sigw1.estimate), \
                random.gauss(params.a2.estimate * self.state[1], params.sigw2.estimate)]
            # below is proportional to the density of two independent gaussians
            # always normalize weights relative to each other, so no need to normalize here
            # removed 1/sqrt(2pi)
            self.weight = self.weight * (params.sigv.estimate**(-1)) * \
             math.exp(-.5 * (self.state[0] - obs[0])**2/params.sigv.estimate**2) * \
             (params.sigv.estimate**(-1)) * \
             math.exp(-.5 * (self.state[1] - obs[1])**2/params.sigv.estimate**2)

        if model == "ar_sigw":
            self.state = random.gauss(self.state, params.sigw.estimate)
            self.weight = self.weight * (2*math.pi*(.5**2))**(-.5) * \
             math.exp(-.5 * (self.state - obs)**2/.5**2)
             # sigv=0.5


        self.state_hist.append(self.state)
        # only need to keep lag+1 elements of state_history
        if len(self.state_hist) == lag+2:
            del self.state_hist[0]



# make particles = rep(particle(),N)

class particles_class:
# contains set of particles, their normalized weights, and ESS

    def __init__(self, model, N, params, observation):
        self.parts = []
        for i in range(N):
            self.parts.append( particle(model, N, params, observation) )
        self.ESS = N

    def update_particle_states(self, model, params, new_obs, lag, sAR_true_params):

        for part in self.parts:
            part.update_particle(model, params, new_obs, lag, sAR_true_params)

    def normalize_particle_weights(self):
        # calculate normalization constant
        sum_weights = 0.0
        for part in self.parts:
            sum_weights += part.weight
        if sum_weights==0:
            print "particle states ", [ part.state for part in self.parts ]
            print "particle weights ", [ part.weight for part in self.parts ]
        # normalize weights, and calculate ESS
        sum_weight_squared = 0.0
        #print "sum_weights:", sum_weights
        for part in self.parts:
            #print "sum_weights:", sum_weights
            part.weight /= sum_weights
            sum_weight_squared += part.weight * part.weight

        self.ESS = 1.0 / sum_weight_squared
        # sanity checking -- normalization constant should now be 1
        sum_weights = 0.0
        for part in self.parts:
            sum_weights += part.weight
        #print "Sum weights now",sum_weights
        ##print "ESS now",self.ESS


    def resample_particles(self):

        num_particles = len(self.parts)
        new_particles = []

        # systematically sample particles
        r = random.uniform(0,1.0/num_particles)
        indx = 0
        w_sum = 0
        # bisect.bisect() would have been easier to get indices
        while r < 1:
            while self.parts[indx].weight + w_sum < r:
                w_sum += self.parts[indx].weight
                indx += 1
            new_particle = copy.deepcopy( self.parts[indx] )
            new_particle.weight = 1.0 / num_particles
            new_particles.append( new_particle )
            r += 1.0/num_particles

        self.parts = new_particles
        self.ESS = num_particles


class param:
# contains information for a single parameter and functions for updating under ioem

    def __init__(self, em_method, in_param):

        self.estimate = in_param
        self.bar = in_param
        if em_method == "ioem":
            self.old_bar = in_param
            self.Sums = Sums()
            self.eta = .5                  # will this cause problems in weighted_regression?
            self.memlen = 1 / self.eta


    # x and y in regression completely different from x and y in particle filters
    def weighted_regression(self, y, full_swing, index, gamma=2, min_memory=500):

        """Function to perform regression"""

        #print "pre Swwwwxx", self.Sums.Swwwwxx
        #print "pre Swwwwx", self.Sums.Swwwwx
        #print "pre Swwww", self.Sums.Swwww
        #print "pre Swwxx", self.Sums.Swwxx
        #print "pre Swwxy", self.Sums.Swwxy
        #print "pre Swwx", self.Sums.Swwx
        #print "pre Swwy", self.Sums.Swwy
        #print "pre Sww", self.Sums.Sww
        #print "memlen", self.memlen

        self.Sums.shift_and_update(self.eta, y)

        #print "Swwx", self.Sums.Swwx
        #print "Swwxx", self.Sums.Swwxx
        #print "Swwwwx", self.Sums.Swwwwx
        #print "Swwwwxx", self.Sums.Swwwwxx
        #print "post Swwwwxx", self.Sums.Swwwwxx
        #print "post Swwwwx", self.Sums.Swwwwx
        #print "post Swwww", self.Sums.Swwww
        #print "post Swwxx", self.Sums.Swwxx
        #print "post Swwxy", self.Sums.Swwxy
        #print "post Swwx", self.Sums.Swwx
        #print "post Swwy", self.Sums.Swwy
        #print "post Sww", self.Sums.Sww

        ### calculate current estimate for the regression

        # normalization constant (1/determinant) of (X'X)^-1
        try:
            norm = 1/(self.Sums.Sww*self.Sums.Swwxx - self.Sums.Swwx*self.Sums.Swwx)
        except:
            #print "Warning: determinant zero"
            norm = 1

        # beta-hat = (X'X)^-1 X'y  (CHECKED)
        beta0 = norm*(self.Sums.Swwxx*self.Sums.Swwy - self.Sums.Swwx*self.Sums.Swwxy)
        beta1 = norm*(-self.Sums.Swwx*self.Sums.Swwy + self.Sums.Sww*self.Sums.Swwxy)

        # s^2 = (X beta - y)'(X beta - y) = beta' X'X beta - 2 y'X beta + y'y  (CHECKED)
        ss = beta0*beta0*self.Sums.Sww + 2*beta0*beta1*self.Sums.Swwx + beta1*beta1*self.Sums.Swwxx \
            -2*(beta0*self.Sums.Swwy + beta1*self.Sums.Swwxy) \
            + self.Sums.Swwyy

        #  var beta-hat = s^2 (X'X)^-1.
        #
        # this formula holds when the weight w_i = sqrt( 1 / variance(y_i) ).
        # Actually, we use weights that do NOT correspond to the variance of each
        # data point, but rather to an a-priori chosen weight sequence.  This means
        # that the variance in beta-hat under the null (each data point has the same
        # variance) is different.  In fact it is:
        #
        # var beta-hat = (X'X)^-1 X' diag(w_i^2 sigma^2) X (X'X)^-1
        #
        # (When w_i == 1 this simplifies to (X'X)^-1).  The matrix X' diag(w_i^2 sigma^2) X
        # has a relatively simple form (list of rows; S[i][j] = S_ij ):
        #
        #   [ [ Swwww, Swwwwx ] , [ Swwwwx, Swwwwxx ] ]

        Xprime_X_inv  = [ [ self.Sums.Swwxx * norm, -self.Sums.Swwx * norm ], [ -self.Sums.Swwx * norm, self.Sums.Sww * norm ] ]
        Xprime_diag_X = [ [ self.Sums.Swwww * ss, self.Sums.Swwwwx * ss ], [ self.Sums.Swwwwx * ss, self.Sums.Swwwwxx * ss ] ]
        varbetahat = matrixProduct( matrixProduct( Xprime_X_inv, Xprime_diag_X ), Xprime_X_inv )

        varbeta0 = varbetahat[0][0]
        varbeta1 = varbetahat[1][1]

        if varbeta0 <= 0:
            #removed as only prints on first run and clogging screen
            #print "Warning: var beta0 = ",varbeta0
            varbeta0 = 1

        if varbeta1 <= 0:
            #removed as only prints on first run and clogging screen
            #print "Warning: var beta1 = ",varbeta1
            varbeta1 = 1


        # calculate memory length, depending on the observed drift
        self.memlen = gamma *  math.sqrt(varbeta0) / ( abs(beta1) + math.sqrt(varbeta1) )

        # Set new weight.  For testing, use uniform weights (1.0 / iteration)
        # Make sure that memory doesn't exceed previous memory
        if full_swing:
            #self.eta = min( index**-0.5 , max(1.0 / self.memlen, index**-1))
            self.eta = min( index**-0.5 , max(1.0 / self.memlen, self.eta / (1.0 + self.eta) ) )
        else:
            # wait for min_mem worth of sufficient statistics before relying on regression
            self.eta = self.eta / (1.0 + self.eta)

        # store eta as memlen - this is recorded for plotting, though not used
        self.memlen = 1.0 / self.eta
            


    def update_est_ioem(self, i, forth, min_mem, sweep_indx, ioem_gam):


        if i>=forth+1+min_mem or sweep_indx>1:
            # full-swing updates
            thetaupdate = (1.0/self.eta) * self.bar + (1 - 1.0/self.eta) * self.old_bar
            self.weighted_regression(thetaupdate, True, index=i, \
             gamma=ioem_gam, min_memory=min_mem)
            
            self.estimate = self.bar    # parameter estimates are updated as in OEM, but with gamma (here eta) set through weighted regression

        else:
            # frozen parameter
            self.weighted_regression(self.bar, False, index=i, \
             gamma=ioem_gam, min_memory=min_mem)




class params:
# contains set of parameters and functions to map from summary sufficient statistics (averaged counts) to estimates

    def __init__(self, model, em_method, in_params):

        if model == "coin":
            self.ff = param(em_method, in_params.ff)
            self.uf = param(em_method, in_params.uf)
            self.ft = param(em_method, in_params.ft)
            self.ut = param(em_method, in_params.ut)

        if model == "AR1":
            self.a = param(em_method, in_params.a)
            self.sigw = param(em_method, in_params.sigw)
            self.sigv = param(em_method, in_params.sigv)

        if model == "sAR1":
            self.sigv = param(em_method, in_params.sigv)

        if model == "sar":
            self.sigv2 = param(em_method, in_params.sigv2)

        if model == "sv":
            self.phi = param(em_method, in_params.phi)
            self.sig = param(em_method, in_params.sig)
            self.b = param(em_method, in_params.b)

        if model == "two_ar":
            self.a1 = param(em_method, in_params.a1)
            self.sigw1 = param(em_method, in_params.sigw1)
            self.sigv1 = param(em_method, in_params.sigv1)
            self.a2 = param(em_method, in_params.a2)
            self.sigw2 = param(em_method, in_params.sigw2)
            self.sigv2 = param(em_method, in_params.sigv2)

        if model == "two_ar_shared_sigv":
            self.a1 = param(em_method, in_params.a1)
            self.sigw1 = param(em_method, in_params.sigw1)
            self.sigv = param(em_method, in_params.sigv)
            self.a2 = param(em_method, in_params.a2)
            self.sigw2 = param(em_method, in_params.sigw2)

        if model == "ar_sigw":
            self.sigw = param(em_method, in_params.sigw)


    # don't go into update params if i<=forth+1, burn-in period
    def update_bars_bem(self, model, counts):

        """Function to update parameters in BEM, only use when i>forth+1"""

        if model=="coin":
            if counts.ff!=0 and counts.f!=0:
                self.ff.bar = counts.ff/counts.f
            if counts.uf!=0 and counts.u!=0:
                self.uf.bar = counts.uf/counts.u
            if counts.ft!=0 and counts.f!=0:
                self.ft.bar = counts.ft/counts.f
            if counts.ut!=0 and counts.u!=0:
                self.ut.bar = counts.ut/counts.u

        if model=="AR1":
            if counts.S1!=0:
                self.a.bar = counts.S2/counts.S1
            # if counts.S1!=0: ???
            #use self.a.bar in defining self.sigw.bar for convenience, eq to S2/S1
            self.sigw.bar = math.sqrt(abs(counts.S3-self.a.bar*counts.S2))
            #if counts.S3-self.a.bar*counts.S2<0:
            #print "S3-S2^2/S1<0, error: NEGATIVE VARIANCE"
            self.sigv.bar = math.sqrt(counts.S4)

        if model=="sAR1":
            self.sigv.bar = math.sqrt(counts.S4)

        if model=="sar":
            self.sigv2.bar = counts.S4

        if model=="sv":
            if counts.S2!=0:
                self.phi.bar = counts.S1/counts.S2
            self.sig.bar = math.sqrt(abs(counts.S3-self.phi.bar*counts.S1))
            #if counts.S3-self.phi.bar*counts.S2<0:
            #print "S3-S1^2/S2<0, error: NEGATIVE VARIANCE"
            self.b.bar = math.sqrt(counts.S4)

        if model=="two_ar":
            if counts.S1_1!=0:
                self.a1.bar = counts.S2_1/counts.S1_1
                self.sigw1.bar = math.sqrt(abs(counts.S3_1-self.a1.bar*counts.S2_1))
                self.sigv1.bar = math.sqrt(counts.S4_1)
            if counts.S1_2!=0:
                self.a2.bar = counts.S2_2/counts.S1_2
                self.sigw2.bar = math.sqrt(abs(counts.S3_2-self.a2.bar*counts.S2_2))
                self.sigv2.bar = math.sqrt(counts.S4_2)

        if model=="two_ar_shared_sigv":
            if counts.S1_1!=0:
                self.a1.bar = counts.S2_1/counts.S1_1
            self.sigw1.bar = math.sqrt(abs(counts.S3_1-self.a1.bar*counts.S2_1))
            if counts.S1_2!=0:
                self.a2.bar = counts.S2_2/counts.S1_2
            self.sigw2.bar = math.sqrt(abs(counts.S3_2-self.a2.bar*counts.S2_2))
            self.sigv.bar = math.sqrt((counts.S4_1+counts.S4_2)/2)

        if model=="ar_sigw":
            self.sigw.bar = math.sqrt(abs(counts.S3-counts.S2**2/counts.S1))

#from counts to theta bars
    # don't go into update bars if i<=forth+1
    def update_bars_oem(self, model, counts):

        """Function to update parameters in OEM, only use when i>forth+1"""

        if model=="coin":
            if counts.ff!=0 and counts.f!=0:
                self.ff.bar = counts.ff/counts.f
            if counts.uf!=0 and counts.u!=0:
                self.uf.bar = counts.uf/counts.u
            if counts.ft!=0 and counts.f!=0:
                self.ft.bar = counts.ft/counts.f
            if counts.ut!=0 and counts.u!=0:
                self.ut.bar = counts.ut/counts.u

        if model=="AR1":
            if counts.S1!=0:
                self.a.bar = counts.S2/counts.S1
            self.sigw.bar = math.sqrt(abs(counts.S3-self.a.bar*counts.S2))
            if counts.S3-self.a.bar*counts.S2<0:
                print "S3-S2^2/S1<0, error: NEGATIVE VARIANCE"
            self.sigv.bar = math.sqrt(counts.S4)

        if model=="sAR1":
            self.sigv.bar = math.sqrt(counts.S4)

        if model=="sar":
            self.sigv2.bar = counts.S4

        if model=="sv":
            if counts.S2!=0:
                self.phi.bar = counts.S1/counts.S2
            self.sig.bar = math.sqrt(abs(counts.S3-self.phi.bar*counts.S1))
            if counts.S3-self.phi.bar*counts.S1<0:
                print "S3-S1^2/S2<0, error: NEGATIVE VARIANCE"
            self.b.bar = math.sqrt(counts.S4)

        if model=="two_ar":
            if counts.S1_1!=0:
                self.a1.bar = counts.S2_1/counts.S1_1
            self.sigw1.bar = math.sqrt(abs(counts.S3_1-self.a1.bar*counts.S2_1))
            self.sigv1.bar = math.sqrt(counts.S4_1)
            if counts.S1_2!=0:
                self.a2.bar = counts.S2_2/counts.S1_2
            self.sigw2.bar = math.sqrt(abs(counts.S3_2-self.a2.bar*counts.S2_2))
            self.sigv2.bar = math.sqrt(counts.S4_2)

        if model=="two_ar_shared_sigv":
            if counts.S1_1!=0:
                self.a1.bar = counts.S2_1/counts.S1_1
            self.sigw1.bar = math.sqrt(abs(counts.S3_1-self.a1.bar*counts.S2_1))
            if counts.S1_2!=0:
                self.a2.bar = counts.S2_2/counts.S1_2
            self.sigw2.bar = math.sqrt(abs(counts.S3_2-self.a2.bar*counts.S2_2))
            self.sigv.bar = math.sqrt((counts.S4_1+counts.S4_2)/2)

        if model=="ar_sigw":
            self.sigw.bar = math.sqrt(abs(counts.S3-counts.S2**2/counts.S1))


    def update_bars_ioem(self, model, counts):
        # need to keep previous theta_bar in order to calculate the update

        if model=="coin":
            self.ff.old_bar = self.ff.bar
            self.uf.old_bar = self.uf.bar
            self.ft.old_bar = self.ft.bar
            self.ut.old_bar = self.ut.bar
            if counts.ff!=0 and counts.f_ff!=0:
                self.ff.bar = counts.ff/counts.f_ff
            if counts.uf!=0 and counts.u_uf!=0:
                self.uf.bar = counts.uf/counts.u_uf
            if counts.ft!=0 and counts.f_ft!=0:
                self.ft.bar = counts.ft/counts.f_ft
            if counts.ut!=0 and counts.u_ut!=0:
                self.ut.bar = counts.ut/counts.u_ut

        if model=="AR1":
            self.a.old_bar = self.a.bar
            self.sigw.old_bar = self.sigw.bar
            self.sigv.old_bar = self.sigv.bar
            if counts.S1a!=0: # keep a.bar the same to avoid /0
                self.a.bar = counts.S2a/counts.S1a
            if counts.S1sigw!=0: # keep a.bar the same to avoid /0
                self.sigw.bar = math.sqrt(abs(counts.S3 - counts.S2sigw**2 / counts.S1sigw))
                if counts.S3 - counts.S2sigw**2 / counts.S1sigw<0:
                    print "S3-S2^2/S1<0 error: NEGATIVE VARIANCE"
            self.sigv.bar = math.sqrt(counts.S4)

        if model == "sAR1":
            self.sigv.old_bar = self.sigv.bar
            self.sigv.bar = math.sqrt(counts.S4)

        if model == "sar":
            self.sigv2.old_bar = self.sigv2.bar
            self.sigv2.bar = counts.S4

        if model=="sv":
            self.phi.old_bar = self.phi.bar
            self.sig.old_bar = self.sig.bar
            self.b.old_bar = self.b.bar
            if counts.S2phi!=0: # keep a.bar the same to avoid /0
                self.phi.bar = counts.S1phi/counts.S2phi
            if counts.S2sig!=0: # keep a.bar the same to avoid /0
                self.sig.bar = math.sqrt(abs(counts.S3 - counts.S1sig**2 / counts.S2sig))
                if counts.S3 - counts.S1sig**2 / counts.S2sig<0:
                    print "S3-S1^2/S2<0 error: NEGATIVE VARIANCE"
            self.b.bar = math.sqrt(counts.S4)

        if model=="two_ar":
            self.a1.old_bar = self.a1.bar
            self.sigw1.old_bar = self.sigw1.bar
            self.sigv1.old_bar = self.sigv1.bar
            if counts.S1a_1!=0: # keep a.bar the same to avoid /0
                self.a1.bar = counts.S2a_1/counts.S1a_1
            if counts.S1sigw_1!=0: # keep a.bar the same to avoid /0
                self.sigw1.bar = math.sqrt(abs(counts.S3_1 - counts.S2sigw_1**2 / counts.S1sigw_1))
            self.sigv1.bar = math.sqrt(counts.S4_1)
            self.a2.old_bar = self.a2.bar
            self.sigw2.old_bar = self.sigw2.bar
            self.sigv2.old_bar = self.sigv2.bar
            if counts.S1a_2!=0: # keep a.bar the same to avoid /0
                self.a2.bar = counts.S2a_2/counts.S1a_2
            if counts.S1sigw_2!=0: # keep a.bar the same to avoid /0
                self.sigw2.bar = math.sqrt(abs(counts.S3_2 - counts.S2sigw_2**2 / counts.S1sigw_2))
            self.sigv2.bar = math.sqrt(counts.S4_2)

        if model=="two_ar_shared_sigv":
            self.a1.old_bar = self.a1.bar
            self.sigw1.old_bar = self.sigw1.bar
            if counts.S1a_1!=0: # keep a.bar the same to avoid /0
                self.a1.bar = counts.S2a_1/counts.S1a_1
            if counts.S1sigw_1!=0: # keep a.bar the same to avoid /0
                self.sigw1.bar = math.sqrt(abs(counts.S3_1 - counts.S2sigw_1**2 / counts.S1sigw_1))
            self.sigv.old_bar = self.sigv.bar
            self.sigv.bar = math.sqrt((counts.S4_1+counts.S4_2)/2)
            self.a2.old_bar = self.a2.bar
            self.sigw2.old_bar = self.sigw2.bar
            if counts.S1a_2!=0: # keep a.bar the same to avoid /0
                self.a2.bar = counts.S2a_2/counts.S1a_2
            if counts.S1sigw_2!=0: # keep a.bar the same to avoid /0
                self.sigw2.bar = math.sqrt(abs(counts.S3_2 - counts.S2sigw_2**2 / counts.S1sigw_2))

        if model=="ar_sigw":
            self.sigw.old_bar = self.sigw.bar
            if counts.S1sigw!=0: # keep a.bar the same to avoid /0
                self.sigw.bar = math.sqrt(abs(counts.S3 - counts.S2sigw**2 / counts.S1sigw))

    #from theta bars to theta ests
    def update_ests_bem(self, model):

        """Function to update parameters in BEM, only use when i>forth+1"""

        if model=="coin":
            self.ff.estimate = self.ff.bar
            self.uf.estimate = self.uf.bar
            self.ft.estimate = self.ft.bar
            self.ut.estimate = self.ut.bar

        if model=="AR1":
            self.a.estimate = self.a.bar
            self.sigw.estimate = self.sigw.bar
            self.sigv.estimate = self.sigv.bar

        if model=="sAR1":
            self.sigv.estimate = self.sigv.bar

        if model=="sar":
            self.sigv2.estimate = self.sigv2.bar

        if model=="sv":
            self.phi.estimate = self.phi.bar
            self.sig.estimate = self.sig.bar
            self.b.estimate = self.b.bar

        if model=="two_ar":
            self.a1.estimate = self.a1.bar
            self.sigw1.estimate = self.sigw1.bar
            self.sigv1.estimate = self.sigv1.bar
            self.a2.estimate = self.a2.bar
            self.sigw2.estimate = self.sigw2.bar
            self.sigv2.estimate = self.sigv2.bar

        if model=="two_ar_shared_sigv":
            self.a1.estimate = self.a1.bar
            self.sigw1.estimate = self.sigw1.bar
            self.a2.estimate = self.a2.bar
            self.sigw2.estimate = self.sigw2.bar
            self.sigv.estimate = self.sigv.bar

        if model=="ar_sigw":
            self.sigw.estimate = self.sigw.bar


    # don't go into update params if i<=forth+1
    def update_ests_oem(self, model):

        """Function to update parameters in OEM, only use when i>forth+1"""

        if model=="coin":
            self.ff.estimate = self.ff.bar
            self.uf.estimate = self.uf.bar
            self.ft.estimate = self.ft.bar
            self.ut.estimate = self.ut.bar

        if model=="AR1":
            self.a.estimate = self.a.bar
            self.sigw.estimate = self.sigw.bar
            self.sigv.estimate = self.sigv.bar

        if model=="sAR1":
            self.sigv.estimate = self.sigv.bar

        if model=="sar":
            self.sigv2.estimate = self.sigv2.bar

        if model=="sv":
            self.phi.estimate = self.phi.bar
            self.sig.estimate = self.sig.bar
            self.b.estimate = self.b.bar

        if model=="two_ar":
            self.a1.estimate = self.a1.bar
            self.sigw1.estimate = self.sigw1.bar
            self.sigv1.estimate = self.sigv1.bar
            self.a2.estimate = self.a2.bar
            self.sigw2.estimate = self.sigw2.bar
            self.sigv2.estimate = self.sigv2.bar

        if model=="two_ar_shared_sigv":
            self.a1.estimate = self.a1.bar
            self.sigw1.estimate = self.sigw1.bar
            self.a2.estimate = self.a2.bar
            self.sigw2.estimate = self.sigw2.bar
            self.sigv.estimate = self.sigv.bar

        if model=="ar_sigw":
            self.sigw.estimate = self.sigw.bar


    def update_ests_ioem(self, model, i, lag, min_mem, sweep_indx, ioem_gam):

        if model=="coin":
            # maybe add parameter constraints here?
            self.ff.update_est_ioem(i, lag, min_mem, sweep_indx, ioem_gam)
            self.uf.update_est_ioem(i, lag, min_mem, sweep_indx, ioem_gam)
            self.ft.update_est_ioem(i, lag, min_mem, sweep_indx, ioem_gam)
            self.ut.update_est_ioem(i, lag, min_mem, sweep_indx, ioem_gam)

        if model=="AR1":
            self.a.update_est_ioem(i, lag, min_mem, sweep_indx, ioem_gam)
            self.sigw.update_est_ioem(i, lag, min_mem, sweep_indx, ioem_gam)
            self.sigv.update_est_ioem(i, lag, min_mem, sweep_indx, ioem_gam)

        if model=="sAR1":
            self.sigv.update_est_ioem(i, lag, min_mem, sweep_indx, ioem_gam)

        if model=="sar":
            self.sigv2.update_est_ioem(i, lag, min_mem, sweep_indx, ioem_gam)

        if model=="sv":
            self.phi.update_est_ioem(i, lag, min_mem, sweep_indx, ioem_gam)
            self.sig.update_est_ioem(i, lag, min_mem, sweep_indx, ioem_gam)
            self.b.update_est_ioem(i, lag, min_mem, sweep_indx, ioem_gam)

        if model=="two_ar":
            self.a1.update_est_ioem(i, lag, min_mem, sweep_indx, ioem_gam)
            self.sigw1.update_est_ioem(i, lag, min_mem, sweep_indx, ioem_gam)
            self.sigv1.update_est_ioem(i, lag, min_mem, sweep_indx, ioem_gam)
            self.a2.update_est_ioem(i, lag, min_mem, sweep_indx, ioem_gam)
            self.sigw2.update_est_ioem(i, lag, min_mem, sweep_indx, ioem_gam)
            self.sigv2.update_est_ioem(i, lag, min_mem, sweep_indx, ioem_gam)

        if model=="two_ar_shared_sigv":
            self.a1.update_est_ioem(i, lag, min_mem, sweep_indx, ioem_gam)
            self.sigw1.update_est_ioem(i, lag, min_mem, sweep_indx, ioem_gam)
            self.a2.update_est_ioem(i, lag, min_mem, sweep_indx, ioem_gam)
            self.sigw2.update_est_ioem(i, lag, min_mem, sweep_indx, ioem_gam)
            self.sigv.update_est_ioem(i, lag, min_mem, sweep_indx, ioem_gam)

        if model=="ar_sigw":
            self.sigw.update_est_ioem(i, lag, min_mem, sweep_indx, ioem_gam)


class Sums:
# contains information for perform weighted regression recursively

    def __init__(self):
        self.Swwwwxx = 0
        self.Swwwwx = 0
        self.Swwww = 1
        self.Swwxx = 0
        self.Swwxy = 0
        self.Swwyy = 0
        self.Swwx = 0
        self.Swwy = 0
        self.Sww = 1

    def shift_and_update(self, weight, y):
        self.Swwwwxx = self.Swwwwxx - 2*self.Swwwwx + self.Swwww
        self.Swwwwx  = self.Swwwwx - self.Swwww
        self.Swwxx   = self.Swwxx - 2*self.Swwx + self.Sww
        self.Swwxy   = self.Swwxy - self.Swwy
        self.Swwx    = self.Swwx - self.Sww

        # Update the counters, using that x' = 0
        wfac = 1.0 - weight
        self.Sww     = wfac * self.Sww + weight    # this is always 1
        self.Swwx    = wfac * self.Swwx
        self.Swwxy   = wfac * self.Swwxy
        self.Swwxx   = wfac * self.Swwxx
        self.Swwy    = wfac * self.Swwy + weight * y
        self.Swwyy   = wfac * self.Swwyy + weight * y * y
        self.Swwww   = wfac*wfac * self.Swwww + weight * weight
        self.Swwwwx  = wfac*wfac * self.Swwwwx
        self.Swwwwxx = wfac*wfac * self.Swwwwxx



class updates:
# contains sufficient statistics

    def __init__(self, model):

        if model == "coin":
            self.ff = 0
            self.uf = 0
            self.ft = 0
            self.ut = 0
            self.f = 0
            self.u = 0

        if model == "AR1":
            self.s1 = 0
            self.s2 = 0
            self.s3 = 0
            self.s4 = 0

        if model == "sAR1" or model == "sar":
            self.s4 = 0

        if model == "sv":
            self.s1 = 0
            self.s2 = 0
            self.s3 = 0
            self.s4 = 0

        if model == "two_ar" or model=="two_ar_shared_sigv":
            #_1 is ar with obs[0],a1,sigw1,sigv1 _2 is ar with obs[1],a2,sigw2,sigv2
            self.s1_1 = 0
            self.s2_1 = 0
            self.s3_1 = 0
            self.s4_1 = 0
            self.s1_2 = 0
            self.s2_2 = 0
            self.s3_2 = 0
            self.s4_2 = 0

        if model == "ar_sigw":
            self.s1 = 0
            self.s2 = 0
            self.s3 = 0


    # only pass in the relevant observation
    def count_new_events(self, model, particles, observation):

        # reset updates to 0
        self.__init__(model)

        if model == "coin":
            #state_hist[0] is X_{t-lag-1}, state_hist[1] is X_{t-lag}
            #part.weight is from time t
            for part in particles.parts:
                if part.state_hist[0] == 'F':
                    self.f += part.weight
                    if part.state_hist[1] == 'F':
                        self.ff += part.weight
                    if observation == 'T':
                        self.ft += part.weight
                else:
                    self.u += part.weight
                    if part.state_hist[1] == 'F':
                        self.uf += part.weight
                    if observation == 'T':
                        self.ut += part.weight

        if model == "AR1":
            for part in particles.parts:
                self.s1 += part.state_hist[0]**2 * part.weight
                self.s2 += part.state_hist[0]*part.state_hist[1] * part.weight
                self.s3 += part.state_hist[1]**2 * part.weight
                self.s4 += (observation - part.state_hist[1])**2 * part.weight

        if model == "sAR1" or model == "sar":
            for part in particles.parts:
                self.s4 += (observation - part.state_hist[1])**2 * part.weight

        if model == "sv":
            for part in particles.parts:
                self.s1 += part.state_hist[0]*part.state_hist[1] * part.weight
                self.s2 += part.state_hist[0]**2 * part.weight
                self.s3 += part.state_hist[1]**2 * part.weight
                self.s4 += (observation**2 * math.exp(-part.state_hist[1])) * part.weight

        if model == "two_ar" or model=="two_ar_shared_sigv":
            for part in particles.parts:
                self.s1_1 += part.state_hist[0][0]**2 * part.weight
                self.s2_1 += part.state_hist[0][0]*part.state_hist[1][0] * part.weight
                self.s3_1 += part.state_hist[1][0]**2 * part.weight
                self.s4_1 += (observation[0] - part.state_hist[1][0])**2 * part.weight
                self.s1_2 += part.state_hist[0][1]**2 * part.weight
                self.s2_2 += part.state_hist[0][1]*part.state_hist[1][1] * part.weight
                self.s3_2 += part.state_hist[1][1]**2 * part.weight
                self.s4_2 += (observation[1] - part.state_hist[1][1])**2 * part.weight

        if model == "ar_sigw":
            for part in particles.parts:
                self.s1 += part.state_hist[0]**2 * part.weight
                self.s2 += part.state_hist[0]*part.state_hist[1] * part.weight
                self.s3 += part.state_hist[1]**2 * part.weight




#double check S1 and S2 aren't mixed up for stochastic volatility model
class counts:
# contains summary sufficient statistics and functions to update them

    # initialized at 0, make sure first update has weight 1
    def __init__(self, model, em_method, in_p):

        if model == "coin":
            self.ff = 0
            self.uf = 0
            self.ft = 0
            self.ut = 0
            if em_method == "ioem":
                self.u_uf = 0
                self.u_ut = 0
                self.f_ff = 0
                self.f_ft = 0
            else:
                self.u = 0
                self.f = 0

        if model == "AR1":
            self.S3 = 0
            self.S4 = 0
            if em_method == "ioem":
                self.S1a = 0
                self.S2a = 0
                self.S1sigw = 0
                self.S2sigw = 0
            else:
                self.S1 = 0
                self.S2 = 0

        if model=="sAR1" or model=="sar":
            self.S4 = 0

        if model=="sv":
            self.S3 = 0
            self.S4 = 0
            if em_method == "ioem":
                self.S1phi = 0
                self.S2phi = 0
                self.S1sig = 0
                self.S2sig = 0
            else:
                self.S1 = 0
                self.S2 = 0

        if model == "two_ar" or model=="two_ar_shared_sigv":
            #_1 first ar, _2 second ar
            self.S3_1 = 0
            self.S4_1 = 0
            if em_method == "ioem":
                self.S1a_1 = 0
                self.S2a_1 = 0
                self.S1sigw_1 = 0
                self.S2sigw_1 = 0
                self.S1a_2 = 0
                self.S2a_2 = 0
                self.S1sigw_2 = 0
                self.S2sigw_2 = 0
            else:
                self.S1_1 = 0
                self.S2_1 = 0
                self.S3_2 = 0
                self.S4_2 = 0
                self.S1_2 = 0
                self.S2_2 = 0

        if model == "ar_sigw":
            self.S3 = 0
            # below is unnecessary for ar_sigw, but keeping for consistency
            if em_method == "ioem":
                self.S1sigw = 0
                self.S2sigw = 0
            else:
                self.S1 = 0
                self.S2 = 0

    #put if i>forth+1 outside of function
    def update_counts_bem(self, model, updates, i ,forth, batch_size):

        j = (i-forth-1) % batch_size + 1

        if model=="coin":
            self.f = self.f * (1-1.0/j) + (1.0/j) * updates.f
            self.u = self.u * (1-1.0/j) + (1.0/j) * updates.u
            self.ff = self.ff * (1-1.0/j) + (1.0/j) * updates.ff
            self.uf = self.uf * (1-1.0/j) + (1.0/j) * updates.uf
            self.ft = self.ft * (1-1.0/j) + (1.0/j) * updates.ft
            self.ut = self.ut * (1-1.0/j) + (1.0/j) * updates.ut

        if model=="AR1":
            self.S1 = self.S1 * (1-1.0/j) + (1.0/j) * updates.s1
            self.S2 = self.S2 * (1-1.0/j) + (1.0/j) * updates.s2
            self.S3 = self.S3 * (1-1.0/j) + (1.0/j) * updates.s3
            self.S4 = self.S4 * (1-1.0/j) + (1.0/j) * updates.s4

        if model=="sAR1" or model=="sar":
            self.S4 = self.S4 * (1-1.0/j) + (1.0/j) * updates.s4

        # could save lines by combining AR1 and sv, but may be more confusing?
        if model=="sv":
            self.S1 = self.S1 * (1-1.0/j) + (1.0/j) * updates.s1
            self.S2 = self.S2 * (1-1.0/j) + (1.0/j) * updates.s2
            self.S3 = self.S3 * (1-1.0/j) + (1.0/j) * updates.s3
            self.S4 = self.S4 * (1-1.0/j) + (1.0/j) * updates.s4

        if model=="two_ar" or model=="two_ar_shared_sigv":
            self.S1_1 = self.S1_1 * (1-1.0/j) + (1.0/j) * updates.s1_1
            self.S2_1 = self.S2_1 * (1-1.0/j) + (1.0/j) * updates.s2_1
            self.S3_1 = self.S3_1 * (1-1.0/j) + (1.0/j) * updates.s3_1
            self.S4_1 = self.S4_1 * (1-1.0/j) + (1.0/j) * updates.s4_1
            self.S1_2 = self.S1_2 * (1-1.0/j) + (1.0/j) * updates.s1_2
            self.S2_2 = self.S2_2 * (1-1.0/j) + (1.0/j) * updates.s2_2
            self.S3_2 = self.S3_2 * (1-1.0/j) + (1.0/j) * updates.s3_2
            self.S4_2 = self.S4_2 * (1-1.0/j) + (1.0/j) * updates.s4_2

        if model=="ar_sigw":
            self.S1 = self.S1 * (1-1.0/j) + (1.0/j) * updates.s1
            self.S2 = self.S2 * (1-1.0/j) + (1.0/j) * updates.s2
            self.S3 = self.S3 * (1-1.0/j) + (1.0/j) * updates.s3


    #put if i>forth+1 outside of function
    def update_counts_oem(self, model, updates, gam):

        if model=="coin":
            self.f = self.f * (1-gam) + (gam) * updates.f
            self.u = self.u * (1-gam) + (gam) * updates.u
            self.ff = self.ff * (1-gam) + (gam) * updates.ff
            self.uf = self.uf * (1-gam) + (gam) * updates.uf
            self.ft = self.ft * (1-gam) + (gam) * updates.ft
            self.ut = self.ut * (1-gam) + (gam) * updates.ut

        if model=="AR1":
            self.S1 = self.S1 * (1-gam) + (gam) * updates.s1
            self.S2 = self.S2 * (1-gam) + (gam) * updates.s2
            self.S3 = self.S3 * (1-gam) + (gam) * updates.s3
            self.S4 = self.S4 * (1-gam) + (gam) * updates.s4

        if model=="sAR1" or model=="sar":
            self.S4 = self.S4 * (1-gam) + (gam) * updates.s4

        if model=="sv":
            self.S1 = self.S1 * (1-gam) + (gam) * updates.s1
            self.S2 = self.S2 * (1-gam) + (gam) * updates.s2
            self.S3 = self.S3 * (1-gam) + (gam) * updates.s3
            self.S4 = self.S4 * (1-gam) + (gam) * updates.s4

        if model=="two_ar" or model=="two_ar_shared_sigv":
            self.S1_1 = self.S1_1 * (1-gam) + (gam) * updates.s1_1
            self.S2_1 = self.S2_1 * (1-gam) + (gam) * updates.s2_1
            self.S3_1 = self.S3_1 * (1-gam) + (gam) * updates.s3_1
            self.S4_1 = self.S4_1 * (1-gam) + (gam) * updates.s4_1
            self.S1_2 = self.S1_2 * (1-gam) + (gam) * updates.s1_2
            self.S2_2 = self.S2_2 * (1-gam) + (gam) * updates.s2_2
            self.S3_2 = self.S3_2 * (1-gam) + (gam) * updates.s3_2
            self.S4_2 = self.S4_2 * (1-gam) + (gam) * updates.s4_2

        if model=="ar_sigw":
            self.S1 = self.S1 * (1-gam) + (gam) * updates.s1
            self.S2 = self.S2 * (1-gam) + (gam) * updates.s2
            self.S3 = self.S3 * (1-gam) + (gam) * updates.s3


    #put if i>forth+1 outside of function
    def update_counts_ioem(self, model, updates, p):

        if model=="coin":
            if p.ff.eta!=0.5:
                self.f_ff = self.f_ff * (1-p.ff.eta) + (p.ff.eta) * updates.f
                self.u_uf = self.u_uf * (1-p.uf.eta) + (p.uf.eta) * updates.u
                self.f_ft = self.f_ft * (1-p.ft.eta) + (p.ft.eta) * updates.f
                self.u_ut = self.u_ut * (1-p.ut.eta) + (p.ut.eta) * updates.u
                self.ff = self.ff * (1-p.ff.eta) + (p.ff.eta) * updates.ff
                self.uf = self.uf * (1-p.uf.eta) + (p.uf.eta) * updates.uf
                self.ft = self.ft * (1-p.ft.eta) + (p.ft.eta) * updates.ft
                self.ut = self.ut * (1-p.ut.eta) + (p.ut.eta) * updates.ut
            else:   #i==lag+1, want to disregard default counts of 0 (use eta=1)
                self.f_ff = updates.f
                self.u_uf = updates.u
                self.f_ft = updates.f
                self.u_ut = updates.u
                self.ff = updates.ff
                self.uf = updates.uf
                self.ft = updates.ft
                self.ut = updates.ut

        if model=="AR1":
            if p.a.eta!=0.5:
                self.S1a = self.S1a * (1-p.a.eta) + (p.a.eta) * updates.s1
                self.S2a = self.S2a * (1-p.a.eta) + (p.a.eta) * updates.s2
                self.S1sigw = self.S1sigw * (1-p.sigw.eta) + (p.sigw.eta) * updates.s1
                self.S2sigw = self.S2sigw * (1-p.sigw.eta) + (p.sigw.eta) * updates.s2
                self.S3 = self.S3 * (1-p.sigw.eta) + (p.sigw.eta) * updates.s3
                self.S4 = self.S4 * (1-p.sigv.eta) + (p.sigv.eta) * updates.s4
            else:   #i==lag+1, want to disregard default counts of 0 (use eta=1)
                self.S1a = updates.s1
                self.S2a = updates.s2
                self.S1sigw = updates.s1
                self.S2sigw = updates.s2
                self.S3 = updates.s3
                self.S4 = updates.s4

        if model=="sAR1":
            if p.sigv.eta!=0.5:
                self.S4 = self.S4 * (1-p.sigv.eta) + (p.sigv.eta) * updates.s4
            else:   #i==lag+1, want to disregard default counts of 0 (use eta=1)
                self.S4 = updates.s4

        if model=="sar":
            if p.sigv2.eta!=0.5:
                self.S4 = self.S4 * (1-p.sigv2.eta) + (p.sigv2.eta) * updates.s4
            else:   #i==lag+1, want to disregard default counts of 0 (use eta=1)
                self.S4 = updates.s4

        if model=="sv":
            if p.phi.eta!=0.5:
                self.S1phi = self.S1phi * (1-p.phi.eta) + (p.phi.eta) * updates.s1
                self.S2phi = self.S2phi * (1-p.phi.eta) + (p.phi.eta) * updates.s2
                self.S1sig = self.S1sig * (1-p.sig.eta) + (p.sig.eta) * updates.s1
                self.S2sig = self.S2sig * (1-p.sig.eta) + (p.sig.eta) * updates.s2
                self.S3 = self.S3 * (1-p.sig.eta) + (p.sig.eta) * updates.s3
                self.S4 = self.S4 * (1-p.b.eta) + (p.b.eta) * updates.s4
            else:   #i==lag+1, want to disreguard default counts of 0 (use eta=1)
                self.S1phi = updates.s1
                self.S2phi = updates.s2
                self.S1sig = updates.s1
                self.S2sig = updates.s2
                self.S3 = updates.s3
                self.S4 = updates.s4

        if model=="two_ar":
            if p.a1.eta!=0.5:
                self.S1a_1 = self.S1a_1 * (1-p.a1.eta) + (p.a1.eta) * updates.s1_1
                self.S2a_1 = self.S2a_1 * (1-p.a1.eta) + (p.a1.eta) * updates.s2_1
                self.S1sigw_1 = self.S1sigw_1 * (1-p.sigw1.eta) + (p.sigw1.eta) * updates.s1_1
                self.S2sigw_1 = self.S2sigw_1 * (1-p.sigw1.eta) + (p.sigw1.eta) * updates.s2_1
                self.S3_1 = self.S3_1 * (1-p.sigw1.eta) + (p.sigw1.eta) * updates.s3_1
                self.S4_1 = self.S4_1 * (1-p.sigv1.eta) + (p.sigv1.eta) * updates.s4_1
            else:   #i==lag+1, want to disregard default counts of 0 (use eta=1)
                self.S1a_1 = updates.s1_1
                self.S2a_1 = updates.s2_1
                self.S1sigw_1 = updates.s1_1
                self.S2sigw_1 = updates.s2_1
                self.S3_1 = updates.s3_1
                self.S4_1 = updates.s4_1
            if p.a2.eta!=0.5:
                self.S1a_2 = self.S1a_2 * (1-p.a2.eta) + (p.a2.eta) * updates.s1_2
                self.S2a_2 = self.S2a_2 * (1-p.a2.eta) + (p.a2.eta) * updates.s2_2
                self.S1sigw_2 = self.S1sigw_2 * (1-p.sigw2.eta) + (p.sigw2.eta) * updates.s1_2
                self.S2sigw_2 = self.S2sigw_2 * (1-p.sigw2.eta) + (p.sigw2.eta) * updates.s2_2
                self.S3_2 = self.S3_2 * (1-p.sigw2.eta) + (p.sigw2.eta) * updates.s3_2
                self.S4_2 = self.S4_2 * (1-p.sigv2.eta) + (p.sigv2.eta) * updates.s4_2
            else:   #i==lag+1, want to disregard default counts of 0 (use eta=1)
                self.S1a_2 = updates.s1_2
                self.S2a_2 = updates.s2_2
                self.S1sigw_2 = updates.s1_2
                self.S2sigw_2 = updates.s2_2
                self.S3_2 = updates.s3_2
                self.S4_2 = updates.s4_2

        if model=="two_ar_shared_sigv":
            if p.a1.eta!=0.5:
                self.S1a_1 = self.S1a_1 * (1-p.a1.eta) + (p.a1.eta) * updates.s1_1
                self.S2a_1 = self.S2a_1 * (1-p.a1.eta) + (p.a1.eta) * updates.s2_1
                self.S1sigw_1 = self.S1sigw_1 * (1-p.sigw1.eta) + (p.sigw1.eta) * updates.s1_1
                self.S2sigw_1 = self.S2sigw_1 * (1-p.sigw1.eta) + (p.sigw1.eta) * updates.s2_1
                self.S3_1 = self.S3_1 * (1-p.sigw1.eta) + (p.sigw1.eta) * updates.s3_1
                self.S4_1 = self.S4_1 * (1-p.sigv.eta) + (p.sigv.eta) * updates.s4_1
            else:   #i==lag+1, want to disregard default counts of 0 (use eta=1)
                self.S1a_1 = updates.s1_1
                self.S2a_1 = updates.s2_1
                self.S1sigw_1 = updates.s1_1
                self.S2sigw_1 = updates.s2_1
                self.S3_1 = updates.s3_1
                self.S4_1 = updates.s4_1
            if p.a2.eta!=0.5:
                self.S1a_2 = self.S1a_2 * (1-p.a2.eta) + (p.a2.eta) * updates.s1_2
                self.S2a_2 = self.S2a_2 * (1-p.a2.eta) + (p.a2.eta) * updates.s2_2
                self.S1sigw_2 = self.S1sigw_2 * (1-p.sigw2.eta) + (p.sigw2.eta) * updates.s1_2
                self.S2sigw_2 = self.S2sigw_2 * (1-p.sigw2.eta) + (p.sigw2.eta) * updates.s2_2
                self.S3_2 = self.S3_2 * (1-p.sigw2.eta) + (p.sigw2.eta) * updates.s3_2
                self.S4_2 = self.S4_2 * (1-p.sigv.eta) + (p.sigv.eta) * updates.s4_2
            else:   #i==lag+1, want to disregard default counts of 0 (use eta=1)
                self.S1a_2 = updates.s1_2
                self.S2a_2 = updates.s2_2
                self.S1sigw_2 = updates.s1_2
                self.S2sigw_2 = updates.s2_2
                self.S3_2 = updates.s3_2
                self.S4_2 = updates.s4_2

        if model=="ar_sigw":
            if p.sigw.eta!=0.5:
                self.S1sigw = self.S1sigw * (1-p.sigw.eta) + (p.sigw.eta) * updates.s1
                self.S2sigw = self.S2sigw * (1-p.sigw.eta) + (p.sigw.eta) * updates.s2
                self.S3 = self.S3 * (1-p.sigw.eta) + (p.sigw.eta) * updates.s3
            else:   #i==lag+1, want to disregard default counts of 0 (use eta=1)
                self.S1sigw = updates.s1
                self.S2sigw = updates.s2
                self.S3 = updates.s3


class records:
# contains information to be stored for analysis

    def __init__(self, model, em_method):

        # could initialize the known eta and memlen for bem and oem

        if model == "coin":
            self.ff_est = []
            self.uf_est = []
            self.ft_est = []
            self.ut_est = []
            if em_method == "ioem":
                self.ff_eta = []
                self.uf_eta = []
                self.ft_eta = []
                self.ut_eta = []
                self.ff_memlen = []
                self.uf_memlen = []
                self.ft_memlen = []
                self.ut_memlen = []

        if model == "AR1":
            self.a_est = []
            self.sigw_est = []
            self.sigv_est = []
            if em_method == "ioem":
                self.a_eta = []
                self.sigw_eta = []
                self.sigv_eta = []
                self.a_memlen = []
                self.sigw_memlen = []
                self.sigv_memlen = []

        if model == "sAR1":
            self.sigv_est = []
            if em_method == "ioem":
                self.sigv_eta = []
                self.sigv_memlen = []

        if model == "sar":
            self.sigv2_est = []
            if em_method == "ioem":
                self.sigv2_eta = []
                self.sigv2_memlen = []

        if model == "sv":
            self.phi_est = []
            self.sig_est = []
            self.b_est = []
            if em_method == "ioem":
                self.phi_eta = []
                self.sig_eta = []
                self.b_eta = []
                self.phi_memlen = []
                self.sig_memlen = []
                self.b_memlen = []

        if model == "two_ar":
            self.a1_est = []
            self.sigw1_est = []
            self.sigv1_est = []
            self.a2_est = []
            self.sigw2_est = []
            self.sigv2_est = []
            if em_method == "ioem":
                self.a1_eta = []
                self.sigw1_eta = []
                self.sigv1_eta = []
                self.a1_memlen = []
                self.sigw1_memlen = []
                self.sigv1_memlen = []
                self.a2_eta = []
                self.sigw2_eta = []
                self.sigv2_eta = []
                self.a2_memlen = []
                self.sigw2_memlen = []
                self.sigv2_memlen = []

        if model == "two_ar_shared_sigv":
            self.a1_est = []
            self.sigw1_est = []
            self.sigv_est = []
            self.a2_est = []
            self.sigw2_est = []
            if em_method == "ioem":
                self.a1_eta = []
                self.sigw1_eta = []
                self.sigv_eta = []
                self.a1_memlen = []
                self.sigw1_memlen = []
                self.sigv_memlen = []
                self.a2_eta = []
                self.sigw2_eta = []
                self.a2_memlen = []
                self.sigw2_memlen = []

        if model == "ar_sigw":
            self.sigw_est = []
            if em_method == "ioem":
                self.sigw_eta = []
                self.sigw_memlen = []

    def add_on(self, model, em_method, par):

        if model == "coin":
            self.ff_est.append(par.ff.estimate)
            self.uf_est.append(par.uf.estimate)
            self.ft_est.append(par.ft.estimate)
            self.ut_est.append(par.ut.estimate)
            if em_method == "ioem":
                self.ff_eta.append(par.ff.eta)
                self.uf_eta.append(par.uf.eta)
                self.ft_eta.append(par.ft.eta)
                self.ut_eta.append(par.ut.eta)
                self.ff_memlen.append(par.ff.memlen)
                self.uf_memlen.append(par.uf.memlen)
                self.ft_memlen.append(par.ft.memlen)
                self.ut_memlen.append(par.ut.memlen)

        if model == "AR1":
            self.a_est.append(par.a.estimate)
            self.sigw_est.append(par.sigw.estimate)
            self.sigv_est.append(par.sigv.estimate)
            if em_method == "ioem":
                self.a_eta.append(par.a.eta)
                self.sigw_eta.append(par.sigw.eta)
                self.sigv_eta.append(par.sigv.eta)
                self.a_memlen.append(par.a.memlen)
                self.sigw_memlen.append(par.sigw.memlen)
                self.sigv_memlen.append(par.sigv.memlen)

        if model == "sAR1":
            self.sigv_est.append(par.sigv.estimate)
            if em_method == "ioem":
                self.sigv_eta.append(par.sigv.eta)
                self.sigv_memlen.append(par.sigv.memlen)

        if model == "sar":
            self.sigv2_est.append(par.sigv2.estimate)
            if em_method == "ioem":
                self.sigv2_eta.append(par.sigv2.eta)
                self.sigv2_memlen.append(par.sigv2.memlen)

        if model == "sv":
            self.phi_est.append(par.phi.estimate)
            self.sig_est.append(par.sig.estimate)
            self.b_est.append(par.b.estimate)
            if em_method == "ioem":
                self.phi_eta.append(par.phi.eta)
                self.sig_eta.append(par.sig.eta)
                self.b_eta.append(par.b.eta)
                self.phi_memlen.append(par.phi.memlen)
                self.sig_memlen.append(par.sig.memlen)
                self.b_memlen.append(par.b.memlen)

        if model == "two_ar":
            self.a1_est.append(par.a1.estimate)
            self.sigw1_est.append(par.sigw1.estimate)
            self.sigv1_est.append(par.sigv1.estimate)
            self.a2_est.append(par.a2.estimate)
            self.sigw2_est.append(par.sigw2.estimate)
            self.sigv2_est.append(par.sigv2.estimate)
            if em_method == "ioem":
                self.a1_eta.append(par.a1.eta)
                self.sigw1_eta.append(par.sigw1.eta)
                self.sigv1_eta.append(par.sigv1.eta)
                self.a1_memlen.append(par.a1.memlen)
                self.sigw1_memlen.append(par.sigw1.memlen)
                self.sigv1_memlen.append(par.sigv1.memlen)
                self.a2_eta.append(par.a2.eta)
                self.sigw2_eta.append(par.sigw2.eta)
                self.sigv2_eta.append(par.sigv2.eta)
                self.a2_memlen.append(par.a2.memlen)
                self.sigw2_memlen.append(par.sigw2.memlen)
                self.sigv2_memlen.append(par.sigv2.memlen)

        if model == "two_ar_shared_sigv":
            self.a1_est.append(par.a1.estimate)
            self.sigw1_est.append(par.sigw1.estimate)
            self.sigv_est.append(par.sigv.estimate)
            self.a2_est.append(par.a2.estimate)
            self.sigw2_est.append(par.sigw2.estimate)
            if em_method == "ioem":
                self.a1_eta.append(par.a1.eta)
                self.sigw1_eta.append(par.sigw1.eta)
                self.sigv_eta.append(par.sigv.eta)
                self.a1_memlen.append(par.a1.memlen)
                self.sigw1_memlen.append(par.sigw1.memlen)
                self.sigv_memlen.append(par.sigv.memlen)
                self.a2_eta.append(par.a2.eta)
                self.sigw2_eta.append(par.sigw2.eta)
                self.a2_memlen.append(par.a2.memlen)
                self.sigw2_memlen.append(par.sigw2.memlen)

        if model == "ar_sigw":
            self.sigw_est.append(par.sigw.estimate)
            if em_method == "ioem":
                self.sigw_eta.append(par.sigw.eta)
                self.sigw_memlen.append(par.sigw.memlen)



### smc functions


def pf(model, obs, N, initial_params, sAR_true_params, em_method, lag=21, batch_size=500, oem_exponent=.75, \
       ioem_gam = 1, sweep_indx = 1, min_mem = 500, single_run = True):
#high level function: performs particle filtering

    ### initialize everything ###
    nSim = len(obs)
    par = params(model, em_method, initial_params)
    new_events = updates(model)
    event_counts = counts(model, em_method, initial_params)
    record = records(model, em_method)
    particles = particles_class(model, N, par, obs[0])


    weight_sum =0
    for part in particles.parts:
        weight_sum += part.weight

    # because particles in sv model are simulated data free, need to normalize their weights:
    if model=="sv":
        particles.normalize_particle_weights()


    return_pars = []
    for i in range(2,nSim+1):

        #print "i={} obs={} a={} mem={} sigw={} mem={} sigv={} mem={}".format(
        #   i,obs[i-1],par.a.estimate,par.a.memlen,par.sigw.estimate,par.sigw.memlen,par.sigv.estimate,par.sigv.memlen)

        if em_method == "online":
            if (i-lag)+(sweep_indx-1)*(nSim-lag-1) > 0:
                gam = 1.0/((i-lag)+(sweep_indx-1)*(nSim-lag-1))**oem_exponent

        ### update particles ###

        particles.update_particle_states(model, par, obs[i-1], lag, sAR_true_params) #obs[i-1] because iteration 2 uses vec element [1]
        particles.normalize_particle_weights()

        ### update counts ###

        if i >= (lag + 1):
            if model == "coin":
                #coin uses the very first observation, other models do not
                new_events.count_new_events(model, particles, obs[i-lag-1])
            else:
                new_events.count_new_events(model, particles, obs[i-lag])

        if i >= (lag + 1):
            if em_method == "batch":
                event_counts.update_counts_bem(model, new_events, i, lag, batch_size)
            elif em_method == "online":
                event_counts.update_counts_oem(model, new_events, gam)
            elif em_method == "ioem":
                event_counts.update_counts_ioem(model, new_events, par)

        ### update theta bars if i>forth+1 ###

        if i >= (lag + 1):
            if em_method == "batch":
                par.update_bars_bem(model, event_counts)
            elif em_method == "online":
                par.update_bars_oem(model, event_counts)
            elif em_method == "ioem":
                par.update_bars_ioem(model, event_counts)

        ### update theta ests if i>forth+1+min_mem ###

        if i >= (lag + 1):                              # events have been recorded
            if em_method == "ioem":
                par.update_ests_ioem(model, i, lag, min_mem, sweep_indx, ioem_gam)
            elif em_method == "online":
                if i >= (lag + 1 + min_mem):            # theta bar has had time to stabilize and can be used for simulations
                    par.update_ests_oem(model)
            elif em_method == "batch":
                if (i - lag) % batch_size == 0:         # we have reached the end of the batch; NEEDS adjustment for multiple sweeps
                    par.update_ests_bem(model)
                    event_counts = counts(model, em_method, initial_params)

        ### resample particles ###

        if particles.ESS < N/2.0:
            particles.resample_particles()


        ### record for plots ###

        ## full record of parameter ests for analyzing specific runs:
        if single_run:
            record.add_on(model, em_method, par)
        else:
            ### record parameters at certain iterations for analyzing lots of runs:
            if i in [500000, 400000, 300000, 200000, 150000, 100000, 90000, 80000, 70000, 60000, 50000, 20000, 10000, 5000, 1000]:
                return_pars.append( copy.deepcopy(par) )

    # Sums are returned with return_pars (par1k.ff.Sums)
    if single_run:
        return record
    else:
        return record_pars


class StartingProb:
# parameters for true X_1 distribution

    def __init__(self, model, p):

        if model == "coin":
            self.f = p

        if model == "AR1":
            self.mu = p[0]
            self.sigma = p[1]

        if model == "sAR1" or model == "sar":
            self.mu = p[0]
            self.sigma = p[1]

        if model == "sv":
            self.mu = p[0]       # 0
            self.var = p[1]    # input.sig**2/(1-input.phi**2)

        if model == "two_ar":
            self.mu1 = p[0]
            self.sigma1 = p[1]
            self.mu2 = p[2]
            self.sigma2 = p[3]

        if model == "two_ar_shared_sigv":
            self.mu1 = p[0]
            self.mu2 = p[1]
            self.sigma = p[2]

        if model == "ar_sigw":
            self.mu = p[0]
            self.sigma = p[1]

class set_params:
# for defining true_parameters, initial_parameters, etc.

    def __init__(self, model, p):
        if model == "coin":
            self.ff = p[0]
            self.uf = p[1]
            self.ft = p[2]
            self.ut = p[3]

        if model == "AR1":
            self.a = p[0]
            self.sigw = p[1]
            self.sigv = p[2]

        if model == "sAR1":
            #not sure why I include p[0] and p[1]
            self.a = p[0]
            self.sigw = p[1]
            self.sigv = p[2]

        if model == "sar":
            #not sure why I include p[0] and p[1]
            self.a = p[0]
            self.sigw = p[1]
            self.sigv2 = p[2]

        if model == "sv":
            self.phi = p[0]
            self.sig = p[1]
            self.b = p[2]

        if model == "two_ar":
            self.a1 = p[0]
            self.sigw1 = p[1]
            self.sigv1 = p[2]
            self.a2 = p[3]
            self.sigw2 = p[4]
            self.sigv2 = p[5]

        if model == "two_ar_shared_sigv":
            self.a1 = p[0]
            self.sigw1 = p[1]
            self.sigv = p[2]
            self.a2 = p[3]
            self.sigw2 = p[4]

        if model == "ar_sigw":
            self.sigw = p[0]

def matrixProduct( mat1, mat2 ):
    # 2x2 matrix product.  mat[i][j] = mat_ij  (ith row, jth column)
    return [ [ mat1[0][0]*mat2[0][0] + mat1[0][1]*mat2[1][0],  mat1[0][0]*mat2[0][1] + mat1[0][1]*mat2[1][1] ],
             [ mat1[1][0]*mat2[0][0] + mat1[1][1]*mat2[1][0],  mat1[1][0]*mat2[0][1] + mat1[1][1]*mat2[1][1] ] ]

