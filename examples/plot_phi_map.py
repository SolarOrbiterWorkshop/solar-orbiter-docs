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

t_start = Time('2022-03-03T09:40', format='isot', scale='utc')
t_end = Time('2022-03-03T9:41', format='isot', scale='utc')

search_results_phi = Fido.search(a.Instrument('PHI'), a.Time(t_start.value, t_end.value), a.soar.Product('phi-hrt-blos'))

"""
To search for two products simultaneously, use the '|' operator:
#results_phi = Fido.search(a.Instrument('PHI'), a.Time(t_start.value, t_end.value), (a.soar.Product('phi-hrt-blos') | a.soar.Product('phi-hrt-icnt')))
"""

print(search_results_phi)

files_phi_hrt = Fido.fetch(search_results_phi[0, 0])#, path='./your/path/to/save/PHI/data/')

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



