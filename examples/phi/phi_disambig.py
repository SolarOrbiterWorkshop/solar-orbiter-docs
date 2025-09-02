r"""
=================================
Disambiguation of SO/PHI data and magnetic field Heliographic components
=================================

This example shows how to use a disambiguation file for SO/PHI magnetic field data. 
The disambigation files (with descriptor ``bamb``) are not yet released in SOAR, but they are all avaible on request. 

After the disambiguation of the magnetic field azimuth, we will compute the Heliographic componenets (`Gary and Hagyard 1990 <https://ui.adsabs.harvard.edu/abs/1990SoPh..126...21G/abstract>`__,  `Sun 2022 <https://ui.adsabs.harvard.edu/abs/2013arXiv1309.2392S/abstract>`__) of the magnetic field, namely :math:`B_{\phi}`, :math:`B_{\theta}`, and :math:`B_r`. 

"""
###############################################
# Download SO/PHI data
# --------------------------------------------
# Dowanload the magnetic field strength (``bmag``), inclination (``binc``), and azimuth (``bazi``) from the SOAR.

import sunpy_soar
from sunpy.net import Fido, attrs as a
from astropy.time import Time
import sunpy.map
import sunpy.visualization.colormaps
import matplotlib.pyplot as plt
import numpy as np
import cmasher as cmr
from astropy.io import fits

import warnings, sunpy
warnings.filterwarnings("ignore", category=sunpy.util.SunpyMetadataWarning)

###############################################
# Loading the data from SOAR

t_start_hrt = Time('2024-03-23T22:29', format='isot', scale='utc')
t_end_hrt = Time('2024-03-23T22:32', format='isot', scale='utc')

search_results_phi_hrt = Fido.search(a.Instrument('PHI'), a.Time(t_start_hrt.value, t_end_hrt.value), a.soar.Product('phi-hrt-bmag') | a.soar.Product('phi-hrt-binc') | a.soar.Product('phi-hrt-bazi'))
print(search_results_phi_hrt)

file_phi_bmag = Fido.fetch(search_results_phi_hrt[0])
file_phi_binc = Fido.fetch(search_results_phi_hrt[1])
file_phi_bazi = Fido.fetch(search_results_phi_hrt[2])

bmag = sunpy.map.Map(file_phi_bmag[0])
binc = sunpy.map.Map(file_phi_binc[0])
bazi = sunpy.map.Map(file_phi_bazi[0])

###############################################
# Disambiguation of the azimuth
# --------------------------------------------
# The ``bamb`` file consists of three maps: 
#  - the first map is a disambiguation bitmap with four methods (ME0 and potential acute, random, or radial potential);
#  - the second map shows where the ME0 method is applied (90 for strong field, 60 for weak field included in the ME0 method by the mask dilation, 50 for weak field where one of the other three methods is applied);
#  - the third map is a confidence map (not implemented yet).

def phi_disambig(bazi,bamb,method=2):
    """
    input
    bazi: magnetic field azimuth. Type: str or array
    bamb: disambiguation fits. Type: str or array
    method: method selected for the disambiguation of not annealed points (0 for potential acute, 1 for random, 2 for radial acute). Type: int (2 as Default)
    
    output
    disbazi: disambiguated azimuth. Type: array
    disambig_mask: mask of the annealed points. Type: array
    """
    if type(bazi) is str:
        bazi = fits.getdata(bazi)
    if type(bamb) is str:
        bamb = fits.getdata(bamb)
    
    flag=bamb[0] 
    
    disambig = flag >> method # x shifted right by n bits as from https://docs.python.org/3/library/stdtypes.html
    disbazi = bazi.copy()
    disbazi[disambig%2 != 0] += 180
    
    disambig_mask = (bamb[1] != 50)

    return disbazi, disambig_mask



bamb_file = "data/solo_L2_phi-hrt-bamb_20240323T223009_V202506041803_0443230201.fits.gz"
disbazi, disambig_mask = phi_disambig(bazi.data, bamb_file)
bazi = sunpy.map.Map((disbazi, bazi.fits_header))

###############################################
# Show maps
# --------------------------------------------

fig, ax = plt.subplots(nrows=1,ncols=3, figsize=(15,5), subplot_kw={"projection":bmag.wcs}, layout='tight')
bmag.plot(axes=ax[0], cmap = 'gnuplot_r', clim=(0,2500), title='Bmag')
binc.plot(axes=ax[1], cmap = cmr.fusion, clim=(0,180), title='Binc')
bazi.plot(axes=ax[2], cmap = 'hsv', clim=(0,360), title='Bazi')

# plt.show()

###############################################
# Heliographic Components
# --------------------------------------------
# Definition of the function to compute the Heliographic components of the magnetic field from the detector frame components.

def phi_b2ptr(bmag, binc, bazi):
    """
    Equivalent Python translation of IDL's phi_b2ptr.pro.
    Args:
        bmag : sunpy.map.Map of the magnetic field stregth (G)
        binc : sunpy.map.Map of the magnetic field inclination (deg)
        bazi : sunpy.map.Map of the magnetic field disambiguated azimuth (deg)
    Returns:
        bptr : 3D numpy array of shape (nx, ny, 3) [Bp, Bt, Br (G);
               Bp is positive when pointing west; Bt is positive when pointing north]
        lonlat : 3D numpy array of shape (nx, ny, 2)
        aa : transformation matrix as a (3,3) array
    """

    # Check dimensions
    ny, nx = bmag.data.shape
    header = bmag.fits_header

    # Convert bvec to B_xi, B_eta, B_zeta
    field = bmag.data
    gamma = np.deg2rad(binc.data)
    psi = np.deg2rad(bazi.data + 90)
    
    b_xi = field * np.sin(gamma) * np.cos(psi)
    b_eta = field * np.sin(gamma) * np.sin(psi)
    b_zeta = field * np.cos(gamma)
    
    # Get Stonyhurst lon/lat
    all_coords = sunpy.map.all_coordinates_from_map(bmag).heliographic_stonyhurst
    phi, lambd = all_coords.lon.value, all_coords.lat.value
    
    lonlat = np.zeros((ny, nx, 2), dtype=np.float32)
    lonlat[:, :, 0] = phi
    lonlat[:, :, 1] = lambd
    
    # Get angles
    b = np.deg2rad(header['CRLT_OBS']) # disk center latitude
    crota = header['CROTA']       
    p = -np.deg2rad(crota) # p-angle
    
    phi0 = header['HGLN_OBS'] # image center longitude
    phi = np.deg2rad(phi - phi0) # phi0 needed to set to 0 the central meridian
    lambd = np.deg2rad(lambd)
    
    sinb = np.sin(b)
    cosb = np.cos(b)
    sinp = np.sin(p)
    cosp = np.cos(p)
    sinphi = np.sin(phi)
    cosphi = np.cos(phi)
    sinlam = np.sin(lambd)
    coslam = np.cos(lambd)
    
    # Compute transformation matrix elements
    aa = np.zeros((ny, nx, 3, 3))
    aa[...,2,0] = coslam * (sinb * sinp * cosphi + cosp * sinphi) - sinlam * cosb * sinp
    aa[...,2,1] = -coslam * (sinb * cosp * cosphi - sinp * sinphi) + sinlam * cosb * cosp
    aa[...,2,2] = coslam * cosb * cosphi + sinlam * sinb
    
    aa[...,1,0] = -sinlam * (sinb * sinp * cosphi + cosp * sinphi) - coslam * cosb * sinp
    aa[...,1,1] = sinlam * (sinb * cosp * cosphi - sinp * sinphi) + coslam * cosb * cosp
    aa[...,1,2] = -sinlam * cosb * cosphi + coslam * sinb
    
    aa[...,0,0] = -sinb * sinp * sinphi + cosp * cosphi
    aa[...,0,1] = sinb * cosp * sinphi + sinp * cosphi
    aa[...,0,2] = -cosb * sinphi
    
    # Apply the transformation
    bptr =  np.einsum('xyij,jxy->ixy', aa, np.array([b_xi, b_eta, b_zeta]))
    
    return bptr, lonlat, aa

bptr, lonlat, aa = phi_b2ptr(bmag, binc, bazi)
bphi = sunpy.map.Map((bptr[0],bmag.fits_header))
btheta = sunpy.map.Map((bptr[1],bmag.fits_header))
br = sunpy.map.Map((bptr[2],bmag.fits_header))

fig, ax = plt.subplots(nrows=1,ncols=3, figsize=(15,5), layout='tight')
ax[0].imshow(bphi.data, cmap = 'hmimag', clim=(-1500,1500), origin='lower'); ax[0].set_title(r'$B_{\phi}$')
ax[1].imshow(btheta.data, cmap = 'hmimag', clim=(-1500,1500), origin='lower');  ax[1].set_title(r'$B_{\theta}$')
ax[2].imshow(br.data, cmap = 'hmimag', clim=(-1500,1500), origin='lower');  ax[2].set_title(r'$B_r$')
for a in ax: a.contour(disambig_mask, color='k', linewidth=0.5)

plt.show()

