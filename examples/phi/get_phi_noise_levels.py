"""
=================================
Get PHI Data Noise Levels
=================================

This example demonstrates how to compute the noise levels in the polarimetric (L2 Stokes) 
and line-of-sight magnetic field (L2 Blos) data from the Solar Orbiter PHI instrument.

"""

import sunpy_soar
from sunpy.net import Fido, attrs as a
from astropy.time import Time
import sunpy.visualization.colormaps
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import scipy.optimize as spo


###############################################################################
# Helper functions
# -----------------------------------------


def find_nearest(array, value):
    """
    return index of nearest value in array to the desired value
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def gaus(x,a,x0,sigma):
    """
    return Gauss function
    """
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def gaussian_fit(a):
    """
    gaussian fit for data 'a' from np.histogram
    """
    xx=a[1][:-1] + (a[1][1]-a[1][0])/2
    y=a[0][:]
    p0=[0.,sum(xx*y)/sum(y),np.sqrt(sum(y * (xx - sum(xx*y)/sum(y))**2) / sum(y))] #weighted avg of bins for avg and sigma inital values
    p0[0]=y[find_nearest(xx,p0[1])-5:find_nearest(xx,p0[1])+5].mean() #find init guess for ampltiude of gauss func
    p,cov=spo.curve_fit(gaus,xx,y,p0=p0)
    return p

def blos_noise(blos_values, show=False):
    """
    plot blos hist + Gaussian fit
    """
    fig = plt.figure(figsize = (8,6))
    hi = plt.hist(blos_values.flatten(), bins=np.linspace(-2e2,2e2,100), histtype='stepfilled', alpha=0.5, edgecolor='black')
    if not show:
        plt.close(fig)
    tmp = [0,0]
    tmp[0] = hi[0].astype('float64')
    tmp[1] = hi[1].astype('float64')

    #guassian fit + label
    p = gaussian_fit(tmp)    
    xx=hi[1][:-1] + (hi[1][1]-hi[1][0])/2
    lbl = r'$\mu=$' + f'{p[1]:.2f} G' + '\n' + r'$\sigma=$' + f'{p[2]:.2f} G'
    if show:
        plt.plot(xx,gaus(xx,*p),'r--', label=lbl)
        plt.xlabel('BLOS [G]')
        plt.ylabel('Counts')
        plt.yscale('log')
        plt.ylim(1e0,np.max(hi[0])*10)
        plt.legend(loc='upper right', fontsize=20)
        plt.tight_layout()
        plt.show()

    return p

def stokes_noise(stokes_values, show=False):
    """
    plot stokes v hist + Gaussian fit
    """
    fig = plt.figure(figsize = (8,6))
    hi = plt.hist(stokes_values.flatten(), bins=np.linspace(-1e-2,1e-2,200), histtype='stepfilled', alpha=0.5, edgecolor='black')
    if not show:
        plt.close(fig)
    tmp = [0,0]
    tmp[0] = hi[0].astype('float64')
    tmp[1] = hi[1].astype('float64')

    #guassian fit + label
    p = gaussian_fit(tmp)    
    
    if show:
        xx=hi[1][:-1] + (hi[1][1]-hi[1][0])/2
        lbl = r'$\mu=$' + f'{p[1]:.5f}'+r' $V/I_c$' + '\n' + r'$\sigma=$' + f'{p[2]:.5f}'+r' $V/I_c$'
        plt.plot(xx,gaus(xx,*p),'r--', label=lbl)
        plt.xlabel(r'Stokes $V/I_c$')
        plt.ylabel('Counts')
        plt.yscale('log')
        plt.ylim(1e0,np.max(hi[0])*10)
        plt.legend(loc='upper right', fontsize=20)
        plt.tight_layout()
        plt.show()

    return p

###############################################################################
# Searching for PHI-HRT Blos and Stokes Data 
# -----------------------------------------
#
# (Everything also applies to FDT data)
# We first search for **Solar Orbiter PHI-HRT** (High Resolution Telescope) **Blos** data
# in a given time range. The search results will return metadata about available files.

t_start_hrt = Time('2024-10-15T18:00', format='isot', scale='utc')
t_end_hrt = Time('2024-10-15T18:05', format='isot', scale='utc')

search_results_phi_hrt = Fido.search(a.Instrument('PHI'), a.Time(t_start_hrt.value, t_end_hrt.value), (a.soar.Product('phi-hrt-stokes') | a.soar.Product('phi-hrt-blos')))
print(search_results_phi_hrt)

sr_stokes = search_results_phi_hrt[0]
sr_blos = search_results_phi_hrt[1]

blos_file = Fido.fetch(sr_blos)#, path='../data/') - testing use case
stokes_file = Fido.fetch(sr_stokes)#, path='../data/')

blos = fits.getdata(blos_file[0])
stokes = fits.getdata(stokes_file[0])

print(blos.shape), print(stokes.shape)

###############################################################################
# Blos noise
# -----------------------------------------
#
# 1. Select an area of Quiet Sun (i.e. no significant magnetic structures)
#     - Avoid the limb (unless you care about that)
#     - Avoid the edges of the FoV if possible (although our example here is alredy cropped so no need to worry)
#
# It will be difficult to avoid some network structures, so there will be remnant signal in the area you choose.
# Hence any 'noise' value we will now generate will be an overestimate.


plt.imshow(blos, cmap='hmimag', vmin=-1500, vmax=1500, origin="lower")
plt.show()

###############################################################################
#
# The region bounded by (x,y) = (0,1200) and (x,y) = (500,1700) is a suitable region for analyzing the Blos noise.
#

x_start, x_end = 0, 500
y_start, y_end = 1200, 1700

region_slice = (slice(y_start, y_end), slice(x_start, x_end))

"""
Here you can see a small peak near -120G, which is most real signal and not noise.
"""

blos_fit = blos_noise(blos[region_slice], show=True)

###############################################################################
# Stokes noise
# -----------------------------------------
#
# Typically when referring to the polarimetric noise, one refers to the noise level in Stokes V at the continuum wavelength.
#
# As before:
#
# 1. Select an area of Quiet Sun (ie no significant magnetic structures)
#     - Avoid the limb (unless you care about that)
#     - Avoid the edges of the FoV
#
# 2. Find and select the continuum wavelength of the Stokes V
#
# To find the continuum wavelength point, look at the WAVEXX keywords in the fits header.
# It is either the 0 or 5 wavelength position
# The continuum wavelength is the one separated by 0.3A from the neighbour file.

hdr = fits.getheader(stokes_file[0])
contpos = hdr['CONTPOS']-1 #header is 1 based index, but python needs 0 based
#WARNING: 'CONTPOS' may not exist in some 2022 or 2023 datasets. They will be reprocessed with updated headers soon.

plt.imshow(stokes[contpos,3,:,:], cmap='gist_heat', vmin=-0.01, vmax = 0.01, origin="lower")
plt.colorbar()
plt.show()

########################################

stokes_fit = stokes_noise(stokes[contpos,3,:,:][region_slice], show=True)