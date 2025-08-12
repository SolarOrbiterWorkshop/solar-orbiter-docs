"""
=================================
Reproject PHI Blos to CEA
=================================

This example demonstrates how to reproject line-of-sight magnetic field
(L2 Blos) data from the Solar Orbiter PHI instrument to a CEA projection.

"""

import sunpy.map
import astropy.units as u
from astropy.wcs import WCS
from astropy.io import fits
import sunpy.coordinates
from sunpy.coordinates import propagate_with_solar_surface
from sunpy.map.header_helper import make_fitswcs_header
import numpy as np
import os
import sunpy_soar
from sunpy.net import Fido, attrs as a
from astropy.time import Time
import sunpy.visualization.colormaps
import matplotlib.pyplot as plt

###############################################################################
# Searching for PHI-HRT Blos Data
# --------------------------------
#
# We first search for **Solar Orbiter PHI-HRT** (High Resolution Telescope) **Blos** data
# in a given time range. The search results will return metadata about available files.

t_start_hrt = Time('2024-10-15T18:00:00', format='isot', scale='utc')
t_end_hrt = Time('2024-10-15T18:00:35', format='isot', scale='utc')

search_results_phi_hrt = Fido.search(a.Instrument('PHI'), a.Time(t_start_hrt.value, t_end_hrt.value), (a.soar.Product('phi-hrt-blos')))
print(search_results_phi_hrt)

sr_phi_hrt = search_results_phi_hrt[0,0]
blos_file = Fido.fetch(sr_phi_hrt)
blos_map = sunpy.map.Map(blos_file[0])

########################################

blos_map.plot_settings['cmap'] = 'hmimag'
blos_map.plot_settings['vmin'] = -1500
blos_map.plot_settings['vmax'] = 1500
blos_map.plot()

###############################################################################
# Reproject PHI Blos to CEA 
# --------------------------------
#
# - Make a CEA WCS header

from sunpy.map.mapbase import SpatialPair
from sunpy.coordinates import frames as cf
from astropy.coordinates import SkyCoord

cea_scale = u.Quantity(SpatialPair(0.03 * u.deg/u.pix, 0.03 * u.deg/u.pix))  # CEA scale in degrees. (same as HMI SHARP maps)

hgs_center = blos_map.center.transform_to(cf.HeliographicStonyhurst(obstime=blos_map.date))

dlon = 0.03 * u.deg # pixel size in longitude
dlat = 0.03 * u.deg # pixel size in latitude
lon_half_width = 25 * u.deg # +/- 30 deg around center of PHI-HRT image
lat_half_height = 15 * u.deg

nx = int(np.round((2 * lon_half_width / dlon).to_value(u.one)))
ny = int(np.round((2 * lat_half_height / dlat).to_value(u.one)))

ref_hgs = SkyCoord(hgs_center.lon, hgs_center.lat, frame=cf.HeliographicStonyhurst, obstime=blos_map.date)

cea_hdr = make_fitswcs_header((ny, nx), ref_hgs, scale=cea_scale, projection_code="CEA",\
                              instrument=blos_map.instrument, wavelength=blos_map.wavelength)

###############################################################################
# Reproject with Hann kernel
# --------------------------------

with propagate_with_solar_surface():
        outmap = blos_map.reproject_to(cea_hdr,algorithm='adaptive', kernel='Hann')

########################################

outmap.plot_settings['cmap'] = 'hmimag'
outmap.plot_settings['vmin'] = -1500
outmap.plot_settings['vmax'] = 1500

fig = plt.figure(figsize=(9,6)) 
ax = plt.subplot(projection=outmap) # WCS-aware axes 
im = outmap.plot(axes=ax)
cbar = fig.colorbar(im, fraction=0.03)
cbar.set_label('BLOS [G]')

###############################################################################
# Remap PHI-HRT BLOS onto HMI SHARP CEA magnetogram 
# --------------------------------
#
# - First get the nearest SHARP dataset
# - Set 'JSOC_EMAIL' environment variable to your registered email address to download data from JSOC
jsoc_email = os.environ["JSOC_EMAIL"]

result = Fido.search(
    a.Time("2024-10-15 18:11:00", "2024-10-15 18:13:00"),
    a.Sample(1*u.hour),
    a.jsoc.Series("hmi.sharp_cea_720s"),
    a.jsoc.PrimeKey("HARPNUM", 12032), #used JSOC, SolarMonitor and HEK to find out the right HARPNUM
    a.jsoc.Notify(jsoc_email),
    a.jsoc.Segment("magnetogram"))
print(result)

########################################

file = Fido.fetch(result)

########################################

sharp_map = sunpy.map.Map(file)
sharp_map.plot_settings['cmap'] = 'hmimag'

fig = plt.figure(figsize=(9,6))
ax = fig.add_subplot(projection=sharp_map)
im = sharp_map.plot(axes=ax)
cbar = fig.colorbar(im, fraction=0.03)
cbar.set_label('BLOS [G]')
plt.show()

###############################################################################
# Plot PHI-HRT CEA Blos (60s) and HMI SHARP CEA (720s) magnetograms side by side
# --------------------------------
#
# - Take a crop and plot.
# - Notice that there is a shift between them. 
# - This can easily be corrected using a cross correlation function such as 'image_register' in the PHI WCS example

submap = outmap.submap(top_right = sharp_map.top_right_coord, bottom_left = sharp_map.bottom_left_coord)
submap.plot_settings = sharp_map.plot_settings #set the same plot settings

########################################

fig = plt.figure(figsize=(12,10))
ax = plt.subplot(211, projection=sharp_map)
im = sharp_map.plot(axes=ax)
cbar = fig.colorbar(im, fraction=0.03)
cbar.set_label('BLOS [G]')

ax1 = plt.subplot(212, projection=sharp_map)
im1 = submap.plot(axes=ax1)
cbar = fig.colorbar(im1, fraction=0.03)
cbar.set_label('BLOS [G]')
plt.tight_layout()
plt.show()

###############################################################################
# Alternatively, reproject the PHI Blos onto the HMI SHARP magnetogram
# --------------------------------

def make_cea_hdr_for_phi(hmi_mag_map, phi_blos_map):
    cea_hdr_phi = make_fitswcs_header(hmi_mag_map.data.shape, \
                                      hmi_mag_map.reference_coordinate.replicate(rsun=phi_blos_map.reference_coordinate.rsun),\
                                      projection_code='CEA',scale=u.Quantity(hmi_mag_map.scale),
                                      instrument = phi_blos_map.instrument,
                                      wavelength = phi_blos_map.wavelength)

    cea_hdr_phi['dsun_obs'] = hmi_mag_map.coordinate_frame.observer.radius.to(u.m).value
    cea_hdr_phi['hglt_obs'] = hmi_mag_map.coordinate_frame.observer.lat.value
    cea_hdr_phi['hgln_obs'] = hmi_mag_map.coordinate_frame.observer.lon.value
    cea_hdr_phi['crpix1'] = hmi_mag_map.fits_header['CRPIX1']
    cea_hdr_phi['crpix2'] = hmi_mag_map.fits_header['CRPIX2']
    cea_hdr_phi['crval1'] = hmi_mag_map.fits_header['CRVAL1']
    cea_hdr_phi['crval2'] = hmi_mag_map.fits_header['CRVAL2']    
    cea_hdr_phi['PC1_1'] = 1
    cea_hdr_phi['PC1_2'] = 0
    cea_hdr_phi['PC2_1'] = 0
    cea_hdr_phi['PC2_2'] = 1
    cea_hdr_phi['cdelt1'] = hmi_mag_map.fits_header['cdelt1']
    cea_hdr_phi['cdelt2'] = hmi_mag_map.fits_header['cdelt2']
    return cea_hdr_phi


def reproject_phi_2_hmi_cea(phi_blos_map, cea_hdr_phi):
    with propagate_with_solar_surface():
        outmap = phi_blos_map.reproject_to(cea_hdr_phi,algorithm='adaptive', kernel='Hann')
    return outmap

########################################

cea_on_hmi_hdr = make_cea_hdr_for_phi(sharp_map, blos_map)
phi_blos_on_sharp_map = reproject_phi_2_hmi_cea(blos_map, cea_on_hmi_hdr)

########################################

fig = plt.figure(figsize=(12,10))
ax = plt.subplot(211, projection=sharp_map)
im = sharp_map.plot(axes=ax)
cbar = fig.colorbar(im, fraction=0.03)
cbar.set_label('BLOS [G]')

phi_blos_on_sharp_map.plot_settings = sharp_map.plot_settings

ax1 = plt.subplot(212, projection=sharp_map)
im1 = phi_blos_on_sharp_map.plot(axes=ax1)
cbar = fig.colorbar(im1, fraction=0.03)
cbar.set_label('BLOS [G]')
plt.tight_layout()
plt.show()