"""
=================================
Finding and Plotting PHI Data
=================================

This example demonstrates how to search for, download, and plot 
Solar Orbiter PHI (Polarimetric and Helioseismic Imager) data using SunPy.

PHI SOAR Product Codes:
----------------------------

**HRT (High Resolution Telescope)**  
 * `phi-hrt-blos` (Line-of-sight Magnetic Field)
 * `phi-hrt-bmag` (Total Magnetic Field Strength)
 * `phi-hrt-binc` (Inclination of Magnetic Field)
 * `phi-hrt-bazi` (Azimuth of Magnetic Field)
 * `phi-hrt-icnt` (Continuum Intensity)
 * `phi-hrt-stokes` (Stokes Parameters)
 * `phi-hrt-vlos` (Line-of-sight Velocity)

**FDT (Full Disc Telescope)**  

Currently, only **Continuum Intensity** (`phi-fdt-icnt`) and **Blos** (`phi-fdt-blos`) are available in SOAR.

Future releases will include:
 * `phi-fdt-bmag`
 * `phi-fdt-binc`
 * `phi-fdt-bazi`
 * `phi-fdt-stokes`
 * `phi-fdt-vlos`
"""

import sunpy_soar
from sunpy.net import Fido, attrs as a
from astropy.time import Time
import sunpy.map
import sunpy.visualization.colormaps
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import matplotlib as mpl


###############################################################################
# Searching for PHI-HRT Blos Data
# --------------------------------
#
# We first search for **Solar Orbiter PHI-HRT** (High Resolution Telescope) **Blos** data
# in a given time range. The search results will return metadata about available files.


t_start_hrt = Time('2024-03-23T20:00', format='isot', scale='utc')
t_end_hrt = Time('2024-03-23T23:59', format='isot', scale='utc')

search_results_phi_hrt = Fido.search(a.Instrument('PHI'), a.Time(t_start_hrt.value, t_end_hrt.value), a.soar.Product('phi-hrt-blos'))

###############################################
# To search for two products simultaneously, use the '|' operator:
# `results_phi = Fido.search(a.Instrument('PHI'), a.Time(t_start.value, t_end.value), (a.soar.Product('phi-hrt-blos') | a.soar.Product('phi-hrt-icnt')))`


print(search_results_phi_hrt)

###############################################################################
# Fetching the First Available PHI-HRT File
# -----------------------------------------
#
# Once we have the search results, we fetch the first available file.
# When a path isnt passed as a kwarg, the files will save locally into sunpy/data. 
# You can also pass `path='./your/path/to/save/PHI/data/')` to choose where to save the data.
files_phi_hrt = Fido.fetch(search_results_phi_hrt[0, 0])


###############################################################################
# Loading and Plotting PHI-HRT Data
# ---------------------------------
#
# The downloaded file is in FITS format. We load it as a `sunpy.map.Map`
# and adjust the plot settings for better visualization.

# Load the downloaded PHI-HRT Blos image
phi_hrt_blos_map = sunpy.map.Map(files_phi_hrt[0])

# Update the Plot settings
phi_hrt_blos_map.plot_settings['cmap'] = 'hmimag'
phi_hrt_blos_map.plot_settings['vmin'] = -1500
phi_hrt_blos_map.plot_settings['vmax'] = 1500

# Plot the PHI image
plt.figure(figsize=(8, 6))
phi_hrt_blos_map.plot()
plt.colorbar()
plt.title("SO/PHI-HRT Blos")
plt.show()

###############################################################################
# Loading PHI-HRT Stokes Data
# -------------------------------------------
#
# The Stokes data product contain the Stokes I,Q,U,V parameters
# They are measued along 5 points in the 6173 Angstrom spectral line, and one in the nearby continuum
# Hence it is a 4D data cube, that cannot as yet be loaded into one sunpy map object.
#
# There is a further detail that the 2022 PHI-HRT Stokes data had the dimensions in a different order to 2023 onwards.
# Below is a function that can visualise both versions at a given wavelength position.

# Get a 2022 example
t_start_hrt = Time('2022-03-07T00:00', format='isot', scale='utc')
t_end_hrt = Time('2022-03-07T00:01', format='isot', scale='utc')

search_results_phi_hrt = Fido.search(a.Instrument('PHI'), a.Time(t_start_hrt.value, t_end_hrt.value), a.soar.Product('phi-hrt-stokes'))
phi_file = Fido.fetch(search_results_phi_hrt[0, 0], path='./')
old_stokes = fits.getdata(phi_file[0])

########################################

# Get a 2023 example
t_start_hrt_new = Time('2023-03-29T11:40', format='isot', scale='utc')
t_end_hrt_new = Time('2023-03-29T11:41', format='isot', scale='utc')

search_results_phi_hrt_new = Fido.search(a.Instrument('PHI'), a.Time(t_start_hrt_new.value, t_end_hrt_new.value), a.soar.Product('phi-hrt-stokes'))
phi_file_new = Fido.fetch(search_results_phi_hrt_new[0, 0], path='./')
new_stokes = fits.getdata(phi_file_new[0])

########################################

# Show the different dimensions
print(f"2022 Stokes shape: {old_stokes.shape} (Y, X, Stokes, Wavelength)")
print(f"2023 Stokes shape: {new_stokes.shape} (Stokes, Wavelength, Y, X)")

###############################################################################
# Helper Functions for Plotting PHI-HRT Stokes Data
# ---------------------------------
#
# Load the HRT field stop map (to mask out the edges)

def load_field_stop(path = None):
    """load hrt field stop

    Parameters
    ----------
    path: str
        location of the field stop file (optional)

    Returns
    -------
    field_stop: numpy ndarray
        the field stop of the HRT telescope.
    """
    if path is None:
        path = "./data/HRT_field_stop.fits"
    
    hdu_list_tmp = fits.open(path)
    field_stop = np.asarray(hdu_list_tmp[0].data, dtype=np.float32)
    field_stop = np.where(field_stop > 0,1,0)
    
    return field_stop

########################################

# Plotting function

def plot_hrt_stokes(stokes,wv,subsec=None, title=None):
    """plot hrt stokes maps at one wavelength

    Parameters
    ----------
    stokes_arr : numpy ndarray
        Full HRT Stokes Array.
    wv : int
        Index for the desired wavelength position.
    subsec: numpy ndarray
        Region of interest to be plotted [start_x,end_x,start_y,end_y]
    title: str
        Title of figure
        
    Returns
    -------
    None
    """
    fs = load_field_stop()[:,::-1]
    fs = np.where(fs > 0,1,np.nan)


    fig, (ax1, ax2) = plt.subplots(2,2, figsize = (15,12))

    assert stokes.ndim == 4
    if stokes.shape[0] == 6 and stokes.shape[1] == 4:
        print('New Stokes order format')
        stokes_arr = stokes
    elif stokes.shape[-2] == 4 and stokes.shape[-1] == 6:
        print('Old Stokes order format')
        stokes_arr = np.transpose(stokes, [3,2,0,1])  # move wavelength to first axis and polarisation to 2nd axis

    fs_start = (2048 - stokes_arr.shape[2])//2
    fs_end = fs_start + stokes_arr.shape[2]
    fs = fs[fs_start:fs_end, fs_start:fs_end]
    
    if subsec is not None:
        start_row, end_row = subsec[2:4]
        start_col, end_col = subsec[:2]
        assert len(subsec) == 4
        assert start_row >= 0 and start_row < 2048
        assert end_row >= 0 and end_row < 2048
        assert start_col >= 0 and start_col < 2048
        assert end_col >= 0 and end_col < 2048
        fs = fs[start_row:end_row, start_col:end_col]
        
    else:
        start_row, start_col = 0,0
        end_row, end_col = stokes_arr.shape[2],stokes_arr.shape[3]

    stokes_arr = stokes_arr[:,:,start_row:end_row,start_col:end_col] * fs[np.newaxis, np.newaxis, :, :]
    
    cmap = mpl.colormaps["gist_heat"]
    cmap.set_bad(color='black')
    
    im1 = ax1[0].imshow(stokes_arr[wv, 0, :,:], cmap = cmap, origin="lower") 
    im2 = ax1[1].imshow(stokes_arr[wv, 1, :,:], cmap = cmap, origin="lower")
    im3 = ax2[0].imshow(stokes_arr[wv, 2, :,:], cmap = cmap, origin="lower") 
    im4 = ax2[1].imshow(stokes_arr[wv, 3, :,:], cmap = cmap, origin="lower")

    fig.colorbar(im1, ax=ax1[0],fraction=0.046, pad=0.04)
    fig.colorbar(im2, ax=ax1[1],fraction=0.046, pad=0.04,ticks=[-0.01, -0.005, 0, 0.005, 0.01])
    fig.colorbar(im3, ax=ax2[0],fraction=0.046, pad=0.04,ticks=[-0.01, -0.005, 0, 0.005, 0.01])
    fig.colorbar(im4, ax=ax2[1],fraction=0.046, pad=0.04,ticks=[-0.01, -0.005, 0, 0.005, 0.01])
    
    clim = 0.01

    im1.set_clim(0, 1.2)
    im2.set_clim(-clim, clim)
    im3.set_clim(-clim, clim)
    im4.set_clim(-clim, clim)

    ax1[0].set_title(r'$I/<I_c>$')
    ax1[1].set_title(r'$Q/<I_c>$')
    ax2[0].set_title(r'$U/<I_c>$')
    ax2[1].set_title(r'$V/<I_c>$')

    ax1[0].text(35,40, '(a)', color = "white", size = 'x-large')
    ax1[1].text(35,40, '(b)', color = "white", size = 'x-large')
    ax2[0].text(35,40, '(c)', color = "white", size = 'x-large')
    ax2[1].text(35,40, '(d)', color = "white", size = 'x-large')
    
    if isinstance(title,str):
         plt.suptitle(title)
    else:
        plt.suptitle(f"SO/PHI-HRT Stokes at Wavelength Index: {wv}")
    plt.tight_layout()
    plt.show()

###############################################################################
# Plot the PHI-HRT Stokes Data
# --------------------------------

# Plot the old stokes

plot_hrt_stokes(old_stokes, wv=3, title='Old Stokes Format') #subsec=[100, 900, 100, 900]

########################################

# Plot the new stokes

plot_hrt_stokes(new_stokes, wv=3, title='New Stokes Format') #subsec=[100, 900, 100, 900]

###############################################################################
# Searching for PHI-FDT Blos Data
# --------------------------------
#
# Now, we repeat the process for **Solar Orbiter PHI-FDT** (Full Disc Telescope) **Blos** data.

t_start_fdt = Time('2024-08-08T00:00', format='isot', scale='utc')
t_end_fdt = Time('2024-08-08T06:00', format='isot', scale='utc')

search_results_phi_fdt = Fido.search(a.Instrument('PHI'), a.Time(t_start_fdt.value, t_end_fdt.value), a.soar.Product('phi-fdt-blos'))
print(search_results_phi_fdt)

files_phi_fdt = Fido.fetch(search_results_phi_fdt[0, 0])#, path='./your/path/to/save/PHI/data/')

###############################################################################
# Loading, Rotating, and Masking PHI-FDT Data
# -------------------------------------------
#
# The PHI-FDT data needs to be **rotated and recentered** for correct visualization.
# We also apply a **mask to remove off-disc pixels**.

phi_fdt_blos_map = sunpy.map.Map(files_phi_fdt[0]).rotate(recenter = True) # Rotate the image to the correct orientation

#clean up the off-disc pixels for better visualization
#here we find the coordinators that are on the solar disk and create a mask
hpc_coords = sunpy.map.all_coordinates_from_map(phi_fdt_blos_map)
mask = ~sunpy.map.coordinate_is_on_solar_disk(hpc_coords)

#create a sunpy map object, with a mask which is applied when plotting
phi_fdt_blos_map = sunpy.map.Map(phi_fdt_blos_map.data,phi_fdt_blos_map.meta, mask=mask)

# Update the Plot settings
phi_fdt_blos_map.plot_settings['cmap'] = 'hmimag'
phi_fdt_blos_map.plot_settings['vmin'] = -1500
phi_fdt_blos_map.plot_settings['vmax'] = 1500

# Plot the PHI image
plt.figure(figsize=(8, 6))
phi_fdt_blos_map.plot()
plt.colorbar()
plt.title("SO/PHI-FDT Blos")
plt.show()
