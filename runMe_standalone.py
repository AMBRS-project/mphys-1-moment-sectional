# This is the main program of the microphysics model

from pylab import *
from get_coag_kernel import get_coag_kernel
from dNdt_all import dNdt_all
#from dNMdt_coag_rev import dNMdt_coag_rev
from scipy.integrate import odeint

################
# Inputs
################

nbins=40
density=1400.
lower_limit=1E-9 #m
upper_limit=1E-5 #m
coag_eff=1.0 # scale_factor for coagulation
temp=293. #K
pres=1.01E5 #pa
num_times=50 # number of output times (for plotting)
tfinal=10.*60.*60. # s

H2SO4=1E8 # Sulfuric acid concentration [cm-3]
NucScale=1E-7 # Scale factor for activation nucleation [s-1]
mfp = 6.5E-8 # mean free path of air [m]
h2so4diff = 1E-5 # diffusivity of h2so4 [m2 s-1]

No = [600.,100.] # initial number in log modes [cm-3]
Dpm = [3E-8,3E-7] # initial median diameter in log modes [m]
sigma = [2.,2.] # initial sigma in log modes

doNucCondCoag = [False,False,True]

################
# Code
################

# initialize bin sizes
dia=logspace(log10(lower_limit),log10(upper_limit),nbins+1)
diam=sqrt(dia[:-1]*dia[1:])
dDp=dia[1:]-dia[:-1]
dlogDp=log10(dia[1:]/dia[:-1])
xk=density*pi/6.*dia**3
xkm=density*pi/6.*diam**3

kna = 2*mfp/diam # knudson #
beta = (1.+kna)/(1.+2.*kna*(1.+kna)) # noncontinuum correction

ftimes=linspace(0.,tfinal,num_times)
btimes=-ftimes

# initialize size distribution
Nk=zeros((nbins))
dNdDp = zeros((nbins)) # number per bin per box
for m in range(0,len(No)):
    dNdDp = dNdDp+No[m]/(sqrt(2.*pi)*(diam)*log(sigma[m])) \
       *exp(-((log(Dpm[m]/diam))**2./ \
       (2.*(log(sigma[m]))**2.)))
Nk=dNdDp*dDp
Nki=Nk
#Mk=Nk*xkm
#Mki=Mk

coag_kernel = get_coag_kernel(temp,pres,density,nbins,\
                              diam,xkm)

inputs=Nk

soln = odeint(dNdt_all,inputs,ftimes,\
         args=(nbins,coag_kernel,coag_eff,xkm,diam,beta,H2SO4,NucScale,doNucCondCoag))


print('Initial total number =', soln[0,:30].sum())
print('Final total number =', soln[-1,:30].sum())
print('Ratio =', soln[-1,:30].sum()/soln[0,:30].sum())

dNdlogDp=soln[:,:]/dlogDp

if dNdlogDp.min() < -1E-6: 
   print('ERROR! non-trival negative values after odeint!')

dNdlogDp[where(dNdlogDp<0.)]=0. # some values were just < 0.

######################
# Make plots
######################
figure()
semilogx(diam*1E6,dNdlogDp[0,:],'b',diam*1E6,dNdlogDp[-1,:],'r',linewidth=1.5)
xlabel('D$_p$ [$\mu$m]')
ylabel('dN/dlogD$_p$ [cm$^{-3}$]')
legend(['Initial','Final'],loc='best')

dNdlogDp=ma.masked_less(dNdlogDp,1E-3)

figure()
pcolor(ftimes,diam*1E6,log10(dNdlogDp).transpose())
ax=gca()
ax.set_yscale('log')
xlabel('Time [s]')
ylabel('D$_p$ [$\mu$m]')
title('dN/dlogD$_p$ [cm$^{-3}$]')
colorbar(format=FormatStrFormatter('$10^{%2.1f}$'))
show()

#pcolor(log10(coag_kernel))
#colorbar()
#show()


#h2so4diff = 1E-5 # diffusivity of h2so4 [m2 s-1]
