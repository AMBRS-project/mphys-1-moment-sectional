from pylab import *
from dNdt_coag import dNdt_coag
from dNdt_cond import dNdt_cond
# This calculates the first order rate of change in particles due to
# different processes.
# units of dNdt returned is particles per cm3 per sec

def dNdt_all(num_in_bin,time,nbins,coag_kernel,coag_eff,bin_mean_mass,diam,
         beta,H2SO4,NucScale,doNucCondCoag):
   
   dNdt=zeros(nbins)

   # Nucleation
   if doNucCondCoag[0]:
      dNdt[0]+=NucScale*H2SO4

   # Condensation
   if doNucCondCoag[1]:
      dNdt+=dNdt_cond(num_in_bin,time,nbins,bin_mean_mass,diam,beta,H2SO4)

   # Coagulation
   if doNucCondCoag[2]:
      dNdt+=dNdt_coag(num_in_bin,time,nbins,coag_kernel,coag_eff,bin_mean_mass)

   return dNdt

