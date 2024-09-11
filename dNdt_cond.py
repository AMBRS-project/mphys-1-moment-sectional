from pylab import *
import numpy

# This calculates the first order rate of change in particles due to
# condensation.
# units of dNdt_net returned is particles per cm3 per sec

def dNdt_cond(num_in_bin,time,nbins,bin_mean_mass,diam,beta,H2SO4):

   h2so4diff = 1E-5 # diffusivity of h2so4 [m2 s-1]
   Av=6.022E23
   M_sulf = 98.079/(Av*1000.)   #kg/molecule, molec weight of sulf acid  #
   Im=2*numpy.pi*diam*h2so4diff*M_sulf*(H2SO4*1E6)*beta # mass growth rate in kg/s/particle
#   print Im/bin_mean_mass
#   print Im
#   print diam
   bmmf=bin_mean_mass+Im # particle mass after 1s of growth
   bmm2=empty((nbins+1))
   bmm2[:-1]=bin_mean_mass
   bmm2[-1]=bmm2[-2]*bmm2[-2]/bmm2[-3]
   dNdt_out=num_in_bin*(bmmf-bmm2[:-1])/(bmm2[1:]-bmm2[:-1])
   dNdt_net=empty((nbins))
   dNdt_net[0]=-dNdt_out[0]
   dNdt_net[1:]=dNdt_out[:-1]-dNdt_out[1:]

   return dNdt_net
