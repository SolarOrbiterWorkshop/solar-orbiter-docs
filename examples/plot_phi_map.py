"""
Finding and Plotting Solar Orbiter PHI Data
===========================================

This example demonstrates how to:
- Search for Solar Orbiter PHI (Polarimetric and Helioseismic Imager) data
- Download a sample image
- Plot it using SunPy

PHI SOAR product Codes:

HRT (High Resolution Telescope)
'phi-hrt-blos'
'phi-hrt-bmag'
'phi-hrt-binc'
'phi-hrt-bazi'
'phi-hrt-icnt'
'phi-hrt-stokes'
'phi-hrt-vlos'

FDT (Full Disc Telescope)

Currently only the continuum Intensity and Line-of-sight Magnetic Field (Blos) have been released to SOAR for FDT:
'phi-fdt-blos'
'phi-fdt-icnt'

In Future additionally:
'phi-fdt-bmag'
'phi-fdt-binc'
'phi-fdt-bazi'
'phi-fdt-stokes'
'phi-fdt-vlos'
"""

import sunpy_soar
from sunpy.net import Fido, attrs as a
from astropy.time import Time
import sunpy.map
import sunpy.visualization.colormaps
import matplotlib.pyplot as plt


# ---------------------------------------------------------------
# HRT Blos map
# ---------------------------------------------------------------

t_start_hrt = Time('2022-03-03T09:40', format='isot', scale='utc')
t_end_hrt = Time('2022-03-03T9:41', format='isot', scale='utc')

search_results_phi_hrt = Fido.search(a.Instrument('PHI'), a.Time(t_start_hrt.value, t_end_hrt.value), a.soar.Product('phi-hrt-blos'))

"""
To search for two products simultaneously, use the '|' operator:
#results_phi = Fido.search(a.Instrument('PHI'), a.Time(t_start.value, t_end.value), (a.soar.Product('phi-hrt-blos') | a.soar.Product('phi-hrt-icnt')))
"""

print(search_results_phi_hrt)

files_phi_hrt = Fido.fetch(search_results_phi_hrt[0, 0])#, path='./your/path/to/save/PHI/data/')

# Load the downloaded PHI-HRT blos image
phi_hrt_blos_map = sunpy.map.Map(files_phi_hrt[0])

# Update the Plot settings
phi_hrt_blos_map.plot_settings['cmap'] = 'hmimag'
phi_hrt_blos_map.plot_settings['vmin'] = -2000
phi_hrt_blos_map.plot_settings['vmax'] = 2000

# Plot the PHI image
plt.figure(figsize=(8, 6))
phi_hrt_blos_map.plot()
plt.colorbar()
plt.title("Solar Orbiter PHI-HRT Blos")
plt.show()

# ---------------------------------------------------------------
# FDT Blos map
# ---------------------------------------------------------------

t_start_fdt = Time('2022-04-08T03:10', format='isot', scale='utc')
t_end_fdt = Time('2022-04-08T03:30', format='isot', scale='utc')

search_results_phi_fdt = Fido.search(a.Instrument('PHI'), a.Time(t_start_fdt.value, t_end_fdt.value), a.soar.Product('phi-fdt-blos'))
print(search_results_phi_fdt)

files_phi_fdt = Fido.fetch(search_results_phi_fdt[0, 0])#, path='./your/path/to/save/PHI/data/')

# Load the downloaded PHI-FDT blos image
phi_fdt_blos_map = sunpy.map.Map(files_phi_fdt[0]).rotate(recenter = True) # Rotate the image to the correct orientation

# Update the Plot settings
phi_fdt_blos_map.plot_settings['cmap'] = 'hmimag'
phi_fdt_blos_map.plot_settings['vmin'] = -2000
phi_fdt_blos_map.plot_settings['vmax'] = 2000

#clean up the off-disc pixels for better visualization

#here we find the coordinators that are on the solar disk and create a mask
hpc_coords = sunpy.map.all_coordinates_from_map(phi_fdt_blos_map)
mask = ~sunpy.map.coordinate_is_on_solar_disk(hpc_coords)

#create a sunpy map object, with a mask which is applied when plotting
phi_fdt_blos_map = sunpy.map.Map(phi_fdt_blos_map.data,phi_fdt_blos_map.meta, mask=mask)

# Plot the PHI image
plt.figure(figsize=(8, 6))
phi_fdt_blos_map.plot()
plt.colorbar()
plt.title("Solar Orbiter PHI-FDT Blos")
plt.show()


