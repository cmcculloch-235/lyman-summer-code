#!/usr/bin/env python
from scipy import *
import Pk_library as PKL
import readlos as RL
import matplotlib.pyplot as plt
import matplotlib.cm as cm

### install pylians on calx with
## python3 -m pip install --user Pylians

## initial conditions setup
pk_camb_ics = loadtxt('data/planck1_matterpower_z=99_rightunits.dat')

## los file name to read in
#los_name = '/data/emergence11/irsic/sherwood/ciaran/out_80_ics/los256_n65536_z99.000.dat'
los_name = '../estimator/out_80_ics/los256_n65536_z99.000.dat'

data = RL.LoSfile(los_name,grid=True,rescale_flux=False,verbose=True)

delta_dm = data.rhoker_DM-1.0
BoxSize = 80.0 ## in Mpc/h
axis = 0
MAS = 'NGP'
threads = 8
verbose = True

# compute power spectrum
Pk_dm = PKL.Pk(delta_dm, BoxSize, axis, MAS, threads, verbose)

# 3D P(k)
k_dm       = Pk_dm.k3D
Pk0_dm     = Pk_dm.Pk[:,0] #monopole
Nmodes_dm  = Pk_dm.Nmodes3D

plt.clf()

print (Pk0_dm*k_dm**3/(2.*pi**2.))

plt.plot(k_dm,Pk0_dm*k_dm**3/(2.*pi**2.),color=cm.tab10.colors[0],label='SpecExtract grid')
plt.plot(10**pk_camb_ics[:,0],10**pk_camb_ics[:,1],color='black',label='ICs')

plt.gca().set_xscale('log')
plt.gca().set_yscale('log')

plt.gca().set_xlim(0.05,1.5e1)

plt.legend().draw_frame(False)
plt.savefig('figures/pk_80_1024_ics.pdf')
