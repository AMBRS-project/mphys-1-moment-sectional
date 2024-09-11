from pylab import *

def get_coag_kernel(temp,pres,density,nbins,bin_mean_dia,bin_mean_mass):

   # This function calculates the coagulation kernel [cm3 s-1]
   
   # convert diameter to meters
   dia = bin_mean_dia
   
   # some values that are constant for the calculations
   gc = 8.314 # Gas constant in J/mol*K
   kb = gc/6.022e23; # boltzman constant
   mu=2.5277e-7*temp**0.75302  # viscosity in kg m-1 s-1
   mfp=2.0*mu/(pres*sqrt(8.0*0.0289/(pi*gc*temp)))  # mean free path of air S&P
                                                    # eqn 8.6 in m
   ck=empty((nbins))
   Dk=empty((nbins))
   for k in range(0,nbins):
      ck[k] = (8*kb*temp/(pi*bin_mean_mass[k]))**(1./2.) # coefficient from s&p table 12.1
      Kn = 2*mfp/dia[k] # Knudson number s&p table 12.1
      Dk[k]=kb*temp/(3.0*pi*mu*dia[k])*((5.0+4.0*Kn+6.0*Kn**2.+\
            18.0*Kn**3.)/(5.0-Kn+(8.0+pi)*Kn**2.));
      # diffusivity table 12.1 s&p in m2/s
   
   # calulate kernel, but don't repeat calculations
   coag_kernel=empty((nbins,nbins))
   for i in range(0,nbins):
      for j in range(0,nbins):
         Kn=4.0*(Dk[i]+Dk[j])/(sqrt(ck[i]**2.+ck[j]**2.)*(dia[i]+dia[j]))  # Knudson
                                    #number for 2 particles S&P eqn 12.51
         beta=(1.0+Kn)/(1.0+2.0*Kn*(1.0+Kn))  #correction factor for
                                              #non-continuum S&P eqn 12.50
         coag_kernel[i,j]=2.0*pi*(dia[i]+dia[j])*(Dk[i]+Dk[j])*beta*1e6 
                    # coagulation coefficient in cm3/s  S&P eqn 12.49

   return coag_kernel
