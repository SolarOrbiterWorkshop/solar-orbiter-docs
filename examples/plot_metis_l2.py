"""
==================================
Finding and Plotting Metis Data
==================================

This example demonstrates how to search for, download, and plot 
Solar Orbiter Metis coronagraph observations using SunPy. This example is adapted from Metis tutorial notebooks avaliable at 
https://github.com/SolarOrbiterWorkshop/solo8_tutorials/tree/main/Metis_tutorial created by Aleksandr Burtovoi et al.

"""

import numpy as np
import matplotlib.pyplot as plt
import sunpy.map
from sunpy.net import Fido, attrs as a
import sunpy_soar
from astropy.io import fits



###############################################################################
# Registering the Solar Orbiter Archive (SOAR)
# --------------------------------------------
#
# By importing `sunpy_soar`, the Solar Orbiter Archive (SOAR) is **automatically**
# registered as a data provider in `Fido`. This allows us to directly query and 
# download Solar Orbiter data just like any other SunPy Fido search.
#
# Now, we will use `Fido` to search for available Metis Level 2 data. Metis L2 data are calibrated to be science-ready and are avaliable in two passbands: 
# visible light (VL) and ultraviolet Lyman-alpha (UV). The VL data are separated into two types: total brightness (tB) and polarized brightness (pB).
#
# More details about Metis data products can be found at http://metis.oato.inaf.it/data_products.html
#
# In this example, we will search for VL total brightness data. The search will return metadata about available files.

# valid search term for Metis products are: 'metis-vl-tb', 'metis-vl-pb', 'metis-uv-image'
search_results = Fido.search(a.Time("2023-04-09 07:30", "2023-04-09 08:30"),
                             a.soar.Product("metis-vl-tb"),
                             a.Level(2))

print(search_results)

###############################################################################
# Fetching the First Available File
# ---------------------------------
#
# The search results contain a list of available files matching our query.
# We can use `Fido.fetch` to download the first available file.

downloaded_files = Fido.fetch(search_results[0, 0])

###############################################################################
# Visualizing the Metis image with `sunpy.map`
# -------------------------------------------
#
# The downloaded file is in **FITS format**, which is commonly used for astronomical
# imaging data. We will use `sunpy.map.Map` to load and display it.
# However, in order to properly visualize the Metis data, we need to modify the metadata of FITS file, cut the inner and outer part of FOV, 
# and mask all the bad pixels.
#
# We will define two functions to do this: `cut_metis_fov` and `read_metis`

def cut_metis_fov(hdu, qflag, mask_value=np.nan):
    """
    Masks regions of the Metis instrument field of view (FOV) in the provided FITS HDU data array.
    This function sets pixels to NaN in the HDU data array based on their location relative to the instrument's
    inner and outer FOV boundaries, as well as a quality flag array. 
    Specifically:
    - Pixels within the inner occulter FOV (dist_iocen < fov1) are masked.
    - Pixels outside the outer FOV (dist_suncen > fov2) are masked.
    - Pixels where the quality flag is zero (qflag == 0) are masked.
    The FOV boundaries are calculated using header information from the HDU, including plate scale and FOV radii.
    Parameters
    ----------
    hdu : astropy.io.fits.PrimaryHDU or similar
        FITS HDU object containing image data and relevant header keywords.
    qflag : numpy.ndarray
        Array of quality flags with the same shape as the HDU data. Pixels with a value of 0 are considered invalid.
    mask_value : float, optional
        Value to use for masking invalid pixels. Default is NaN. 
        Use mask_value = 0 to mask with zeros, which works better for several image enhancement algorithms.
    Raises
    ------
    ValueError
        If the plate scale in the x and y directions (CDELT1 and CDELT2) are not equal.
    Modifies
    --------
    hdu.data : numpy.ndarray
        The data array in the HDU is modified in-place, with masked regions set to NaN.
    """

    # check if plate scale is the same in x and y direction
    if hdu.header['CDELT1'] != hdu.header['CDELT2']:
        raise ValueError("Error. CDELT1 != CDELT2 for {fname}".format(fname=hdu.header['FILENAME']))
    # Get FOV in pixel
    fov1 = hdu.header['INN_FOV']*3600/hdu.header['CDELT1']  # pix
    fov2 = hdu.header['OUT_FOV']*3600/hdu.header['CDELT2']  # pix
    # Create meshgrid of pixel coordinates
    x = np.arange(0, hdu.header['NAXIS1'], 1)
    y = np.arange(0, hdu.header['NAXIS2'], 1)
    xx, yy = np.meshgrid(x, y, sparse=True)
    # Calculate distance from Sun center and occulter center
    dist_suncen = np.sqrt((xx-hdu.header['SUN_XCEN'])**2 + (yy-hdu.header['SUN_YCEN'])**2)
    dist_iocen = np.sqrt((xx-hdu.header['IO_XCEN'])**2 + (yy-hdu.header['IO_YCEN'])**2)
    # Mask data outside FOV and bad pixels
    hdu.data[dist_iocen < fov1] = mask_value
    hdu.data[dist_suncen > fov2] = mask_value
    # Mask bad pixels based on quality flag
    hdu.data[qflag == 0] = mask_value

def read_metis(filepath, rot=True):
    """
    Reads a Solar Orbiter Metis FITS file and returns a SunPy map object.
    This function opens the FITS file, processes the data to cut the field of view (FOV)
    and mark bad pixels, and then creates a SunPy map from the processed data.
    Parameters
    ----------
    filepath : str
        Path to the Metis FITS file.
    rot : bool, optional
        If True, the resulting SunPy map will be rotated to have solar north up. Default is True.
    Returns         
    -------
    sunpy.map.Map
        A SunPy map object containing the processed Metis data.
    """
    hdu0 = fits.open(filepath)[0]
    hdu1 = fits.open(filepath)[1]
    cut_metis_fov(hdu0, hdu1.data)
    hdu0.header['RSUN_OBS'] = hdu0.header['RSUN_ARC']
    map_metis = sunpy.map.Map(hdu0.data, hdu0.header)
    if rot == True:
        map_metis = map_metis.rotate()

    return map_metis

###############################################################################
# Now we can use the `read_metis` function to read the downloaded Metis FITS file
# and create a SunPy map. We can then plot the map using its built-in plot method.

metis_vl_tb_file = downloaded_files[0]
metis_vl_tb_map = read_metis(metis_vl_tb_file)

# Define a colormap that handles NaN values (bad pixels) - can be changed to any colormap
Metis_VL_CMAP = plt.get_cmap('afmhot').copy()
Metis_VL_CMAP.set_bad(color='tab:gray')  # np.nan values are in gray

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(projection=metis_vl_tb_map)
im = metis_vl_tb_map.plot(axes=ax, cmap=Metis_VL_CMAP)
metis_vl_tb_map.draw_limb(axes=ax, color='white', linewidth=1.0)
fig.colorbar(im, label='Mean Solar Brightness (MSB)')
ax.set_title('Metis VL tB '+ metis_vl_tb_map.date.strftime('%Y-%m-%d %H:%M:%S'))
plt.show()

###############################################################################
# Note that the time of observation is s/c time, which is different from Earth time.
# To get Earth time, access the 'date_ear' dict in the map metadata.



