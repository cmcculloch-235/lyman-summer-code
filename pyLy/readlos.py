#!/usr/bin/env python

import numpy as np

def rreplace(s, old, new, count):
    return (s[::-1].replace(old[::-1], new[::-1], count))[::-1]

class LoSfile(object):
    def __init__(self,los_name,grid=False,rescale_flux=False,obs=2,verbose=False):
        self.los_name = los_name
        self.tau_name = rreplace(los_name,'los','tau',1)
        self.get_flux_field(verbose=verbose)
        if (verbose):
            self.print_info()
        self.reshape_arrays(grid=grid,verbose=verbose)
        self.rescale_tau(rescale_flux=rescale_flux,obs=obs,verbose=verbose)
        
    def get_flux_field(self,los_name=None,verbose=False):
        if (los_name!=None):
            self.los_name = los_name
            self.tau_name = rreplace(los_name,'los','tau',1)

        if (verbose):
            print ("LoS file: ",self.los_name)
        with open(self.los_name,'rb') as file_ptr:
            self.ztime = np.fromfile(file_ptr,dtype=np.float64,count=1)[0]  ## Redshift
            self.omegam = np.fromfile(file_ptr,dtype=np.float64,count=1)[0] ## Omega_m = Omegda_c + Omega_b
            self.omegal = np.fromfile(file_ptr,dtype=np.float64,count=1)[0]     ## Omega_Lambda
            self.omegab = np.fromfile(file_ptr,dtype=np.float64,count=1)[0]     ## Omega_baryon
            self.h100 = np.fromfile(file_ptr,dtype=np.float64,count=1)[0]       ## h = H0/(100 km/s/Mpc)
            self.box100 = np.fromfile(file_ptr,dtype=np.float64,count=1)[0]     ## Lbox [ckpc/h]
            self.Xh = np.fromfile(file_ptr,dtype=np.float64,count=1)[0]         ## Hydrogen mass fraction
            self.nbins = np.fromfile(file_ptr,dtype=np.int32,count=1)[0]        ## number of bins along LoS
            self.nlos = np.fromfile(file_ptr,dtype=np.int32,count=1)[0]         ## number of LoS
    
            ## E.g. if LoS is in x direction, then ilos=1 and x=0.5*Lbox; 
            ## (y,z) tells you the position in the plane perpendicular to the LoS
            ## LoS axis direction (1=x,2=y,3=z)
            self.ilos = np.fromfile(file_ptr,dtype=np.int32,count=self.nlos)         
            self.xlos = np.fromfile(file_ptr,dtype=np.float64,count=self.nlos)       ## LoS x position
            self.ylos = np.fromfile(file_ptr,dtype=np.float64,count=self.nlos)       ## LoS y position
            self.zlos = np.fromfile(file_ptr,dtype=np.float64,count=self.nlos)       ## LoS z position

            self.posaxis = np.fromfile(file_ptr,dtype=np.float64,count=self.nbins)   ## position along LoS [ckpc/h]
            self.velaxis = np.fromfile(file_ptr,dtype=np.float64,count=self.nbins)   ## position along LoS [km/s]

            ## gas overdensity, Delta=rho/rho_crit
            self.rhoker_H = np.fromfile(file_ptr,dtype=np.float64,count=self.nbins*self.nlos) 
            ## neutral hydrogen fraction, n_HI/n_H
            self.rhoker_H1 = np.fromfile(file_ptr,dtype=np.float64,count=self.nbins*self.nlos) 
            ## gas temperature T [K], HI weighted
            self.tempker_H1 = np.fromfile(file_ptr,dtype=np.float64,count=self.nbins*self.nlos) 
            ## gas velocity v_pec [km/s], HI weighted
            self.velker_H1 = np.fromfile(file_ptr,dtype=np.float64,count=self.nbins*self.nlos) 

            ## DM overdensity, Delta=rho/rho_crit=np/ntot
            self.rhoker_DM = np.fromfile(file_ptr,dtype=np.float64,count=self.nbins*self.nlos) 
            ## gas overdensity, Delta=rho/rho_crit=np/ntot
            self.rhoker_gas = np.fromfile(file_ptr,dtype=np.float64,count=self.nbins*self.nlos) 

        if (verbose):
            print ("tau file: ",self.tau_name)
        with open(self.tau_name,'rb') as file_ptr:
            ## HI optical depth
            self.tau_H1 = np.fromfile(file_ptr,dtype=np.float64,count=self.nbins*self.nlos)

    def reshape_arrays(self,grid=False,verbose=False):
        ## by default reshape into nlos x nbins
        ## if grid==True reshape into sqrt(nlos) x sqrt(nlos) x nbins and assume float32
        if (grid==True):
            n1 = int(np.sqrt(self.nlos))
            if (verbose):
                print ("Reshaping onto 3d grid: ",n1,n1,self.nbins)
            self.rhoker_H = np.float32(self.rhoker_H)
            self.rhoker_H1 = np.float32(self.rhoker_H1)
            self.tempker_H1 = np.float32(self.tempker_H1)
            self.velker_H1 = np.float32(self.velker_H1)
            self.rhoker_DM = np.float32(self.rhoker_DM)
            self.rhoker_gas = np.float32(self.rhoker_gas)
            self.tau_H1 = np.float32(self.tau_H1)

            self.rhoker_H = np.reshape(self.rhoker_H,((n1,n1,self.nbins)))
            self.rhoker_H1 = np.reshape(self.rhoker_H1,((n1,n1,self.nbins)))
            self.tempker_H1 = np.reshape(self.tempker_H1,((n1,n1,self.nbins)))
            self.velker_H1 = np.reshape(self.velker_H1,((n1,n1,self.nbins)))
            self.rhoker_DM = np.reshape(self.rhoker_DM,((n1,n1,self.nbins)))
            self.rhoker_gas = np.reshape(self.rhoker_gas,((n1,n1,self.nbins)))
            self.tau_H1 = np.reshape(self.tau_H1,((n1,n1,self.nbins)))
        else:
            if (verbose):
                print ("Reshaping onto los grid: ",self.nlos,self.nbins)
            self.rhoker_H = np.reshape(self.rhoker_H,((self.nlos,self.nbins)))
            self.rhoker_H1 = np.reshape(self.rhoker_H1,((self.nlos,self.nbins)))
            self.tempker_H1 = np.reshape(self.tempker_H1,((self.nlos,self.nbins)))
            self.velker_H1 = np.reshape(self.velker_H1,((self.nlos,self.nbins)))
            self.rhoker_DM = np.reshape(self.rhoker_DM,((self.nlos,self.nbins)))
            self.rhoker_gas = np.reshape(self.rhoker_gas,((self.nlos,self.nbins)))
            self.tau_H1 = np.reshape(self.tau_H1,((self.nlos,self.nbins)))

    def print_info(self):
        print ("Redshift,Omega_m,Omega_l,Omega_b:")    
        print ("  ",self.ztime,self.omegam,self.omegal,self.omegab)
    
        Hz = 100.0*self.h100*np.sqrt(self.omegam*(1.+self.ztime)**3. + self.omegal)
        print ("Hubble rate H(z) [km/s/Mpc]:",Hz)
    
        vmax = self.box100*1e-3 * Hz/self.h100 * 1./(1.+self.ztime)
        dvbin = vmax/self.nbins
        print ("Boxsize [ckpc/h]:",self.box100)    
        print ("Boxsize [km/s]",vmax)
        print ("Bin size [km/s]: ",dvbin)

    def rescale_tau(self,rescale_flux=True,obs=2,verbose=False):
        ##############################################
        ## rescale mean transmission
        ##############################################
        
        if (obs==0):
            ## Kim et al. 2005, High-res Lya P1d
            taueff = 0.0023*(1.+self.ztime)**3.65
        elif (obs==1):
            ## Palanque-Delabroiulle et al. 2013, BOSS Lya P1d
            taueff = 0.0046*(1.+self.ztime)**3.3
        elif (obs==2):
            ## from Viel+13, high-res Lya P1d + SDSS-II (Becker et al. 2011)
            if (self.ztime < 4.4):
                taueff = -0.132 + 0.751 * ((1.0 + self.ztime)/(1.0 + 3.5))**2.90 # below 4.4
            else:
                taueff = 1.142 * ((1.0 + self.ztime)/(1.0 + 4.4))**4.91 ## above and including 4.4 
        elif (obs==3):
            ## Boera+19, high-res Lya P1d (z>4)
            taueff = 1.56*((1.+self.ztime)/5.75)**4.

        flux_obs = np.exp(-taueff) * self.nlos*self.nbins
        count=0
        tol=1e-8
        scale=1.0
        newscale=2.0

        if (rescale_flux==False):
            if (verbose):
                print ('Redshift = ',self.ztime)
                print ('Original mean H1 flux = ',np.mean(np.exp(-self.tau_H1)))
                print ('Observed mean H1 flux = ',np.exp(-taueff))
            self.flux = np.exp(-self.tau_H1)
            return
    
        while (fabs(newscale-scale)>tol*newscale):
            count=count+1
            scale=newscale
            mean_flux=0
            tau_mean_flux=0

            temp=np.exp(-scale*self.tau_H1)
            mean_flux=sum(temp)
            tau_mean_flux=sum(temp*self.tau_H1)

            newscale = scale + (mean_flux-flux_obs)/tau_mean_flux
            if (newscale < 0): newscale = 0

        if (verbose):
            print ('Redshift = ',self.ztime)
            print ('Original mean H1 flux = ',np.mean(np.exp(-self.tau_H1)))
        self.tau_H1 = self.tau_H1 * scale
        if (verbose):
            print ('Rescaled mean H1 flux = ',np.mean(np.exp(-self.tau_H1)))
            print ('Observed mean H1 flux = ',np.exp(-taueff))
            print ('scale = ',scale)
        self.flux = np.exp(-scale*self.tau_H1)
