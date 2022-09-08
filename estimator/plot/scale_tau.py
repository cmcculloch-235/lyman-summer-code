#!/usr/bin/env python

from scipy import *
import numpy as np

data_dir="../out_80_snap18_new/"
data_z = 2

los_name = data_dir + "/los256_n65536_z" + str(data_z)  + ".000.dat"
tau_name = data_dir + "/tau256_n65536_z" + str(data_z)  + ".000.dat"

with open(los_name,'rb') as file_ptr:
    ztime = np.fromfile(file_ptr,dtype=float64,count=1)[0]      ## Redshift
    omegam = np.fromfile(file_ptr,dtype=float64,count=1)[0]     ## Omega_m = Omegda_c + Omega_b
    omegal = np.fromfile(file_ptr,dtype=float64,count=1)[0]     ## Omega_Lambda
    omegab = np.fromfile(file_ptr,dtype=float64,count=1)[0]     ## Omega_baryon
    h100 = np.fromfile(file_ptr,dtype=float64,count=1)[0]       ## h = H0/(100 km/s/Mpc)
    box100 = np.fromfile(file_ptr,dtype=float64,count=1)[0]     ## Lbox [ckpc/h]
    Xh = np.fromfile(file_ptr,dtype=float64,count=1)[0]         ## Hydrogen mass fraction
    nbins = np.fromfile(file_ptr,dtype=int32,count=1)[0]        ## number of bins along LoS
    nlos = np.fromfile(file_ptr,dtype=int32,count=1)[0]         ## number of LoS
    
    ## E.g. if LoS is in x direction, then ilos=1 and x=0.5*Lbox; 
    ##      (y,z) tells you the position in the plane perpendicular to the LoS
    ilos = np.fromfile(file_ptr,dtype=int32,count=nlos)         ## LoS axis direction (1=x,2=y,3=z)
    xlos = np.fromfile(file_ptr,dtype=float64,count=nlos)       ## LoS x position
    ylos = np.fromfile(file_ptr,dtype=float64,count=nlos)       ## LoS y position
    zlos = np.fromfile(file_ptr,dtype=float64,count=nlos)       ## LoS z position

    posaxis = np.fromfile(file_ptr,dtype=float64,count=nbins)   ## position along LoS [ckpc/h]
    velaxis = np.fromfile(file_ptr,dtype=float64,count=nbins)   ## position along LoS [km/s]

    ## gas overdensity, Delta=rho/rho_crit
    rhoker_H = np.fromfile(file_ptr,dtype=float64,count=nbins*nlos) 
    ## neutral hydrogen fraction, n_HI/n_H
    rhoker_H1 = np.fromfile(file_ptr,dtype=float64,count=nbins*nlos) 
    ## gas temperature T [K], HI weighted
    tempker_H1 = np.fromfile(file_ptr,dtype=float64,count=nbins*nlos) 
    ## gas velocity v_pec [km/s], HI weighted
    velker_H1 = np.fromfile(file_ptr,dtype=float64,count=nbins*nlos) 

with open(tau_name,'rb') as file_ptr:
    ## HI optical depth
    tau_H1 = np.fromfile(file_ptr,dtype=float64,count=nbins*nlos)

####################
tau_H1 = np.reshape(tau_H1,((nlos,nbins)))
####################


"""
print((tau_H1[0]))
exit(1)
####################
### Choose LoS
jlos = 0 ## number in [0,nlos-1]
print(("IGM properties along LoS ",jlos))
directions = ['x','y','z']
print(("  direction of LoS:",directions[ilos[jlos]-1]))
pos_los = [xlos[jlos],ylos[jlos],zlos[jlos]]; del pos_los[ilos[jlos]-1]
print(("  position of LoS [ckpc/h]: ",pos_los[0],pos_los[1]))

## tau field
tau_arr = tau_H1[jlos,:] * 1.0
print(("Opitcal depth along LoS ",jlos))
print(("  ",tau_arr))
print(("  min/max: ",min(tau_arr),max(tau_arr)))
with open('python/Output_0_Optical_Depth_Sims.txt','w') as file_ptr:
    for x in tau_arr:
        s = '{0:.16e}'.format(x)
        file_ptr.write(s+'\n')
"""

def rescale_UVB(redshift,tau_in,UVB_option="Viel",output_flux=True,rescale=True):
    ##############################################
    ########### Rescale the np.mean flux ############
    ##############################################

    z = redshift
    if (UVB_option == "Kim"):
        print('Using Kim+05 observations!')
        taueff = 0.0023*(1.+redshift)**3.65 # Kim+05
    elif (UVB_option == "P-D"):
        taueff = 0.0046*(1.+z)**3.3 # P-D et al BOSS Lya 1D 2013
    elif (UVB_option == "Boera"):
        taueff = 1.56*((1.+z)/5.75)**4. # Boera et al. high res Lya 1D 2019
    elif (UVB_option == "Viel"):
        # Viel et al 2013 Lya 1D                                                                    
        if (z <= 4.5):
            taueff = 0.751*((1.+z)/4.5)**2.90 - 0.132
        else:
            taueff = 2.26*((1.+z)/6.2)**4.91
        
    flux_lya_H1_obs = np.exp(-taueff)

    # Rescale using Newton-Raphson iteration
    if (rescale):
        #default start: 0.01
        scale = 1e-1
        i = 0
        while i < 1000: ## this might be better written with while statement
            i += 1
            scale_old = scale
            print(i, scale, end=" ")
            scale = scale + (np.mean(np.exp(-scale*tau_in)) - flux_lya_H1_obs)/np.mean(tau_in*np.exp(-scale*tau_in))
            dscale = abs(scale - scale_old)
            print(dscale)
            if (dscale < 1e-8):
                break

        print('Redshift = ',redshift)
        print('Original np.mean H1 flux = ',np.mean(np.exp(-tau_in)))
        tau_out = tau_in * scale
        print('Rescaled np.mean H1 flux = ',np.mean(np.exp(-tau_out)))
        print('Observed np.mean H1 flux = ',flux_lya_H1_obs)
        print('scale = ',scale)
    else:
        print('Redshift = ',redshift)
        print('Original np.mean H1 flux = ',np.mean(np.exp(-tau_in)))
        tau_out = tau_in * 1.0
    if (output_flux==True):
        # compute the transmitted flux
        return np.exp(-tau_out)
    else:
        return tau_out




tau_out = rescale_UVB(ztime, tau_H1, UVB_option="Viel", output_flux=False)
print("write...")

tau_out.tofile(data_dir + "/tau_r.dat")
