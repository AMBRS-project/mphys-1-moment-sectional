from pylab import *

# This calculates the first order rate of change in particles due to
# coagulation.
# units of coag_ch returned is particles per cm3 per sec

def dNdt_coag(num_in_bin,time,nbins,coag_kernel,coag_eff,bin_mean_mass):

   coag_ch = zeros((nbins))
   for i in range(0,nbins):
      for j in range(i,nbins):
         if i == j:
            scale = 1./2.
         else:
            scale = 1.
         coag_rate = scale*coag_eff*coag_kernel[i,j]*num_in_bin[i]*num_in_bin[j] 
         # the number of coagulation events between i and j per second per cm3
         coag_ch[i] = coag_ch[i]-coag_rate
         coag_ch[j] = coag_ch[j]-coag_rate
         coag_mass=bin_mean_mass[i]+bin_mean_mass[j]
         # find which bin this goes in
         k=j
         while ((k < nbins-1) and (bin_mean_mass[k+1] < coag_mass)):
            k=k+1
         if (k < nbins-2):
            fc=(coag_mass-bin_mean_mass[k])/(bin_mean_mass[k+1]-bin_mean_mass[k])
            coag_ch[k]=coag_ch[k]+coag_rate*(1-fc) # add one particle per cm3
                   # from bin per coagulation event
            coag_ch[k+1]=coag_ch[k+1]+coag_rate*fc
         elif (k==nbins-2):
            fc=(coag_mass-bin_mean_mass[k])/(bin_mean_mass[k+1]-bin_mean_mass[k])
            coag_ch[k]=coag_ch[k]+coag_rate*(1-fc)

   return coag_ch
