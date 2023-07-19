import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import os
import math

from scipy.stats import norm, gamma

'''
sns.set_context("poster")
sns.set_style("ticks")

mpl.rc('text', usetex = True)
mpl.rc('font', family = 'serif')

'''

class Mcmc_3D:

    def __init__(self, bounds, displacements, prior_rvs, filename, phys_model):

        # Possible values of random variable X
        self.xbounds  = np.array(bounds[0])
        self.ybounds  = np.array(bounds[1])
        self.zbounds  = np.array(bounds[2])
        
        # Initialize a random number generator
        self.prng     = np.random.default_rng()

        # Max displacements for every degree of freedom
        self.displacements = displacements

        # Store the list of prior random variables as an attribute
        self.prior_rvs = prior_rvs

        # Load measured data from file
        self.ind_var, self.measurement = np.loadtxt(filename)

        # Function object that helps compute what reponse should be for a
        # particular value of the independent variable
        self.physical_model = phys_model

        
    def get_log_prior(self, u):

        # Compute sum of log priors
        log_prior = 0
        
        for i, p_rv in enumerate(self.prior_rvs):
            log_prior += math.log(p_rv.pdf(u[i]))
        
        
        return log_prior


    def get_log_likelihood(self, u):

        # Computes log likelihood from data
        
        nmeas           = self.ind_var.size
        log_likelihood  = 0
        
        for i in range(nmeas):

            # Mean should be whatever the physical model predicts
            mean = self.physical_model(self.ind_var[i], u)

            # Likelihood of observing data i
            log_likelihood += math.log(norm(loc = mean, scale = u[2]).pdf(self.measurement[i]))

        return log_likelihood

    
    def get_log_posterior(self, u):

        # Get log posterior from log prior and log likelihood

        return self.get_log_prior(u) + self.get_log_likelihood(u)

    
    # Get value of log f at state u where f is the pdf
    def get_log_f_val(self, u):
        
        # Compute log of posterior
        return self.get_log_posterior(u)
    
    # Random walk
    def get_next_pos (self, u_cur):
        
        u_new = np.copy(u_cur)

        # Pick one of the three random variables
        i_var = self.prng.choice(3)

        # Displace this random variable
        u_new[i_var] = u_cur[i_var] + self.displacements[i_var]*self.prng.uniform (low = -1, high = 1)
        
        return u_new
        

    def check_bounds(self, u_cur):

        # Returns False if walker is out of bounds
        
        if u_cur[0] > self.xbounds[-1] or u_cur[0] < self.xbounds[0]:
            return False
        
        if u_cur[1] > self.ybounds[-1] or u_cur[1] < self.ybounds[0]:
            return False

        if u_cur[2] > self.zbounds[-1] or u_cur[2] < self.zbounds[0]:
            return False

        return True


    def make_plot(self, states):

        # Makes a plot of the trajectory
        
        fig, ax1 = plt.subplots(figsize = (7,5))

        ax2 = ax1.twinx()
        
        n, a = states.shape
        
        ax1.plot(10*np.arange(n)+1, states[:,0])
        ax2.plot(10*np.arange(n)+1, states[:,1])
        ax2.plot(10*np.arange(n)+1, states[:,2])

        ax1.set_ylabel(r'$\sigma$')
        ax2.set_ylabel(r'$a,b$')
        ax1.set_xlabel('MCS')

        plt.tight_layout()
        plt.savefig('trajec.pdf')
        plt.show()
        
    def run_mc(self, u0, nsteps = 200):

        # Main MC loop
        u_old          = np.array(u0)
        
        states = []
        
        for istep in range(nsteps):

            # Compute f at current u
            log_f_old = self.get_log_f_val(u_old)

            # Propose a new u
            u_new  = self.get_next_pos (u_old)

            # Reject the move if bounds are exceeded
            if not self.check_bounds(u_new):
                continue
            
            # Compute f at new u
            log_f_new = self.get_log_f_val (u_new)

            del_log_f = log_f_new - log_f_old
            
            # Compute acceptance factor
            # if del_log_f is larger than zero accenptance factor will be one

            if del_log_f > 0:
                acc_factor = 1
            else:
                acc_factor = math.exp(log_f_new - log_f_old)

            # Proposal probabilities are equal both ways

            # Accept/Reject proposal
            if (self.prng.random() <= acc_factor):

                # Update position
                u_old = u_new
            
            if ((istep+1) % 10 == 0):
                states.append(u_old)
                
                print('MC Step # : ',istep + 1, ' Parameter: ', u_old)
            print(u_old)
        # Uncomment this line out if you want to make a plot 
        self.make_plot(np.array(states))
                
        # Return estimates

        
        return 

if __name__=='__main__':

    # u is an array u[0] is a and u[1] is b
    def activation ( voltage, u ):
        #change this code for the physical model using A and inverse Kappa instead
        act = 1 / (1 + math.exp((voltage - u[0])*(-u[1])))
        return act

    # Priors on three parameters a, b and sigma
    # Define three random variables    
    # Change them to the A and 1/k
    rv_a = norm(0, scale = 100)
    rv_b = norm(0, scale = 2)
    rv_s = gamma(6.25, scale = 1/125)

    # List of prior random variables
    priors = [rv_a, rv_b, rv_s]
    
    # Range of each parameter
    bounds = [[-100,100],[-5, 5],[0,1]]

    # Maximum displacement for each parameter
    disps = [5, 0.05, 0.02]
    
    # Create a Mcmc_3D object
    mcmc_obj = Mcmc_3D(bounds, disps, priors, 'expt_obs_I.dat', activation)

    # Run length
    nsteps   = 1000

    # Run mcmc trajectory
    mcmc_obj.run_mc([10, 0, 0.05], nsteps)
    
    
    
