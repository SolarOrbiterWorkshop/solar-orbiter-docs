"""
=================================
WCS coordinates in SO/PHI-HRT Data
=================================

This example shows how to check the WCS in the SO/PHI-HRT header and how to possibly correct them. 
The data released to SOAR since 2025 hasbetter WCS thanks to a correlation with SO/PHI-FDT data, when possible.

The source of the error in the WCS is the focus mechanism of SO/PHI-HRT, the boresight of the instrument, and the small variations fo the latter with temperature and pointing, which is not fully understood yet.

This example shows how to correct the WCS using magnetograms, but it can be run with continuum images too. 
Other quantities, such as the magnetic field strength, could work, but some functions would need to be re-adapted. 
"""
###############################################
# Download SO/PHI-HRT magnetogram
# --------------------------------------------
# The easiest way to check the WCS of SO/PHI-HRT is by overlaying the data on a SO/PHI-FDT map or, if available, on a SDO/HMI map (also on JHelioviewer).

import sunpy_soar
from sunpy.net import Fido, attrs as a
from astropy.time import Time
import sunpy.map
import sunpy.visualization.colormaps
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from sunpy.coordinates import SphericalScreen
import datetime
import os

import warnings, sunpy
warnings.filterwarnings("ignore", category=sunpy.util.SunpyMetadataWarning)

###############################################
# Loading the data from SOAR

t_start_hrt = Time('2024-03-23T20:00', format='isot', scale='utc')
t_end_hrt = Time('2024-03-23T23:59', format='isot', scale='utc')

search_results_phi_hrt = Fido.search(a.Instrument('PHI'), a.Time(t_start_hrt.value, t_end_hrt.value), a.soar.Product('phi-hrt-blos'))

files_phi_hrt = Fido.fetch(search_results_phi_hrt[0, 0])
hrt_map = sunpy.map.Map(files_phi_hrt[0])
original_header = hrt_map.fits_header.copy()

###############################################
# Correction of CROTA by 0.15 degrees (known value)

def rotate_header(h,angle,center = [1024.5,1024.5]):
    """calculate new header when image is rotated by a fixed angle

    Parameters
    ----------
    h : astropy.io.fits.header.Header
        header of image to be rotated
    angle : float
        angle to rotate image by (in degrees)
    center : list or numpy.array
        center of rotation (x,y) in pixel, Default is [1024.5,1024.5]

    Returns
    -------
    h: astropy.io.fits.header.Header
        new header
    """
    if 'CROTA' in h:
        k = 'CROTA' # positive angle = clock-wise rotation of the reference system axes 
    else:
        k = 'CROTA2'
    h[k] -= angle
    h['PC1_1'] = np.cos(h[k]*np.pi/180)
    h['PC1_2'] = -np.sin(h[k]*np.pi/180)
    h['PC2_1'] = np.sin(h[k]*np.pi/180)
    h['PC2_2'] = np.cos(h[k]*np.pi/180)
    rad = angle * np.pi/180
    rot = np.asarray([[np.cos(rad),-np.sin(rad),0],[np.sin(rad),np.cos(rad),0],[0,0,1]])
    coords = np.asarray([h['CRPIX1'],h['CRPIX2'],1])
#     center = [1024.5,1024.5] # CRPIX from 1 to 2048, so 1024.5 is the center
    tr = np.asarray([[1,0,center[0]],[0,1,center[1]],[0,0,1]])
    invtr = np.asarray([[1,0,-center[0]],[0,1,-center[1]],[0,0,1]])

    M = tr @ rot @ invtr
    bl = M @ np.asarray([0,0,1])
    tl = M @ np.asarray([0,2048,1])
    br = M @ np.asarray([2048,0,1])
    tr = M @ np.asarray([2048,2048,1])

    O = -np.asarray([bl,tl,br,tr]).min(axis=0)[:-1]
    newO = np.asarray([[1,0,O[0]+1],[0,1,O[1]+1],[0,0,1]])
    newM = newO @ M
    new_coords = newM @ coords
    h['CRPIX1'] = round(new_coords[0],4)
    h['CRPIX2'] = round(new_coords[1],4)
    
    return h

def center_coord(hdr):
    """calculate the center of the solar disk in the rotated reference system

    Parameters
    ----------
    hdr : header
        header of the fits file

    Returns
    -------
    center: [x,y,1] coordinates of the solar disk center (units: pixel)
    """
    pxsc = hdr['CDELT1']
    crval1 = hdr['CRVAL1']
    crval2 = hdr['CRVAL2']
    crpix1 = hdr['CRPIX1']
    crpix2 = hdr['CRPIX2']
    if 'PC1_1' in hdr:
        PC1_1 = hdr['PC1_1']
        PC1_2 = hdr['PC1_2']
        PC2_1 = hdr['PC2_1']
        PC2_2 = hdr['PC2_2']
    else:
        if 'CROTA2' in hdr:
            CROTA = hdr['CROTA2']
        else:
            CROTA = hdr['CROTA']
        PC1_1 = np.cos(CROTA*np.pi/180)
        PC1_2 = -np.sin(CROTA*np.pi/180)
        PC2_1 = np.sin(CROTA*np.pi/180) 
        PC2_2 = np.cos(CROTA*np.pi/180)
    
    HPC1 = 0
    HPC2 = 0
    
    x0 = crpix1 + 1/pxsc * (PC1_1*(HPC1-crval1) - PC1_2*(HPC2-crval2)) - 1
    y0 = crpix2 + 1/pxsc * (PC2_2*(HPC2-crval2) - PC2_1*(HPC1-crval1)) - 1
    
    return np.asarray([x0,y0,1])

crota_manual_correction=0.15
h_hrt = rotate_header(original_header.copy(),-crota_manual_correction, center=center_coord(original_header))
hrt_map = sunpy.map.Map((hrt_map.data,h_hrt))

###############################################
# Download the closest SDO/HMI magnetogram in time
# --------------------------------------------
# Replace ``os.environ["JSOC_EMAIL"]`` with your email address registered on JSOC, or set it in your .bashrc file as `export JSOC_EMAIL=your.account@email.com`.

def downloadClosestHMI(ht,jsoc_email,path=False,cad='45'):
    """
    Script to download the HMI m_45 or ic_45 cosest in time to the provided SO/PHI observation.
    TAI convention and light travel time are taken into consideration.
    
    Parameters
    ----------
    ht: astropy.io.fits.header.Header
        header of the SO/PHI observation
    jsoc_email: str
        email address to be used for JSOC connection
    path: bool
        if True, the path of the cache directory and of the HMI dataset will return as output (DEFAULT: False)
    cad: str
        cadence of the HMI data to be downloaded ('45' and '720' are accepted, DEFAULT: '45')

    Returns
    -------
    hmi_map: sunpy.ma.Map
        HMI map of the closest data set to the SO/PHI observation
    cache_dir: str
        path to sunpy cache directory (if path = True)
    hmi_name: str
        path to the downloaded SDO/HMI file (if path = True)
    """
    
    import drms
    import sunpy, sunpy.map
    from astropy.constants import c
    from astropy import units as u

    t_obs = datetime.datetime.fromisoformat(ht['DATE-AVG'])
    dtai = datetime.timedelta(seconds=37) # datetime.timedelta(seconds=94)
    
    if type(cad) != str:
        cad = str(int(cad))
    if cad == '45':
        dcad = datetime.timedelta(seconds=35) # half HMI cadence (23) + margin
    elif cad == '720':
        dcad = datetime.timedelta(seconds=360+60) # half HMI cadence (23) + margin
    else:
        print('wrong HMI cadence, only 45 and 720 are accepted')
        return None
    
    dltt = datetime.timedelta(seconds=ht['EAR_TDEL']) # difference in light travel time S/C-Earth

    kwlist = ['T_REC','T_OBS','DATE-OBS','CADENCE','DSUN_OBS']
    
    client = drms.Client(email=jsoc_email) 

    lt = np.nan
    n = 0
    
    while np.isnan(lt):
        n += 2
        # not implemented for 720s cadence specific data products, yet!
        if ht['PHIDTYPE'] == 'blos':
            keys = client.query('hmi.m_'+cad+'s['+(t_obs+dtai-dcad+dltt).strftime('%Y.%m.%d_%H:%M:%S')+'-'+
                            (t_obs+dtai+dcad+dltt).strftime('%Y.%m.%d_%H:%M:%S')+']',seg=None,key=kwlist,n=n)
        elif ht['PHIDTYPE'] == 'vlos':
            keys = client.query('hmi.v_'+cad+'s['+(t_obs+dtai-dcad+dltt).strftime('%Y.%m.%d_%H:%M:%S')+'-'+
                            (t_obs+dtai+dcad+dltt).strftime('%Y.%m.%d_%H:%M:%S')+']',seg=None,key=kwlist,n=n)
        else:
            keys = client.query('hmi.ic_'+cad+'s['+(t_obs+dtai-dcad+dltt).strftime('%Y.%m.%d_%H:%M:%S')+'-'+
                            (t_obs+dtai+dcad+dltt).strftime('%Y.%m.%d_%H:%M:%S')+']',seg=None,key=kwlist,n=n)
        keys = keys[keys['T_OBS'] != 'MISSING']
        if np.size(keys['T_OBS']) > 0:
            lt = (np.nanmean(keys['DSUN_OBS'])*u.m - ht['DSUN_OBS']*u.m)/c
        else:
            print('adding 60s margin')
            dcad += datetime.timedelta(seconds=60)
        
    dltt = datetime.timedelta(seconds=lt.value) # difference in light travel time S/C-SDO

    T_OBS = [(ind,np.abs((datetime.datetime.strptime(t,'%Y.%m.%d_%H:%M:%S_TAI') - dtai - dltt - t_obs).total_seconds())) for ind, t in zip(keys.index,keys['T_OBS'])]
    ind = T_OBS[np.argmin([t[1] for t in T_OBS])][0]

    if ht['PHIDTYPE'] == 'blos':
        name_h = 'hmi.m_'+cad+'s['+keys['T_REC'][ind]+']{Magnetogram}'
    elif ht['PHIDTYPE'] == 'vlos':
        name_h = 'hmi.v_'+cad+'s['+keys['T_REC'][ind]+']{Dopplergram}'
    else:
        name_h = 'hmi.ic_'+cad+'s['+keys['T_REC'][ind]+']{Continuum}'

    if np.abs((datetime.datetime.strptime(keys['T_OBS'][ind],'%Y.%m.%d_%H:%M:%S_TAI') - dtai - dltt - t_obs).total_seconds()) > np.ceil(int(cad)/2):
        print('WARNING: Closer file exists but has not been found.')
        print(name_h)
        print('T_OBS:',datetime.datetime.strptime(keys['T_OBS'][ind],'%Y.%m.%d_%H:%M:%S_TAI') - dtai - dltt)
        print('DATE-AVG:',t_obs)
        print('')
    else:
        print('HMI T_OBS (corrected for TAI and Light travel time):',datetime.datetime.strptime(keys['T_OBS'][ind],'%Y.%m.%d_%H:%M:%S_TAI') - dtai - dltt)
        print('PHI DATE-AVG:',t_obs)
    
    s45 = client.export(name_h,protocol='fits')
    hmi_map = sunpy.map.Map(s45.urls.url[0],cache=False)
    cache_dir = sunpy.data.CACHE_DIR+'/'
    hmi_name = cache_dir + s45.urls.url[0].split("/")[-1]

    if path:
        return hmi_map, cache_dir, hmi_name
    else:
        return hmi_map

hmi_map = downloadClosestHMI(hrt_map.fits_header,os.environ["JSOC_EMAIL"],path=False,cad='45')

###############################################
# Plotting the results of the Downloads
# --------------------------------------------
plt.figure(figsize=(11,6))
ax = plt.subplot(121,projection=hmi_map)
hmi_map.plot(clim=(-1500,1500),cmap='hmimag', axes=ax)
hmi_map.draw_limb(color='r')
hrt_map.draw_limb(color='w')
constant_lon = SkyCoord(hrt_map.observer_coordinate.lon, np.linspace(-90, 90, 20) * u.deg,
                        frame=hmi_map.observer_coordinate)
plt.gca().plot_coord(constant_lon, color="yellow")
with SphericalScreen(hrt_map.observer_coordinate, only_off_disk=True):
    hrt_map.draw_extent(axes=ax, color='k', lw=1, ls='dotted')
hrt_map.draw_extent(axes=ax, color='k', lw=1)

ax = plt.subplot(122,projection=hrt_map)
hrt_map.plot(clim=(-1500,1500),cmap='hmimag', axes=ax)

plt.tight_layout()
plt.show()

###############################################
# Overplot HMI on SO/PHI-HRT
# --------------------------------------------
from matplotlib.widgets import Slider
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Create the figure and axes
fig, ax = plt.subplots(figsize=(8,8))
plt.subplots_adjust(top=0.85)

# Display the first image
img1 = ax.imshow(hrt_map.data, cmap='hmimag', clim=(-1500,1500), origin='lower')
hmi_remap = hmi_map.reproject_to(hrt_map.wcs)
# Display the second image with initial alpha
img2 = ax.imshow(hmi_remap.data, cmap='hmimag', clim=(-1500,1500), origin='lower', alpha=0.5)

# Create slider axis and slider
slider_ax = inset_axes(ax, width="60%", height="5%", loc='upper center', borderpad=0.6)
slider = Slider(slider_ax, '', 0.0, 1.0, valinit=0.5)
label_text = slider_ax.text(0.5, 1.6, 'Transparency [HRT --> HMI]', transform=slider_ax.transAxes,
                            ha='center', va='bottom', fontsize=10)

# Update function for the slider
def update(val):
    alpha = slider.val
    img2.set_alpha(alpha)
    fig.canvas.draw_idle()

slider.on_changed(update)

plt.show()

###############################################
# The maps are clearly not aligned, you can see two negative (blue) sunspots in the center of the image instead of one.

###############################################
# Check if the WCS by HRT have been updated
# --------------------------------------------
# HISTORY will show if the WCS have been updated or not

hist = hrt_map.fits_header['HISTORY']

for line in hist:
    if 'WCS' in line:
        print(line)

###############################################
# Functions needed
# --------------------------------------------
# We define a remap function to keep the header infomration we want when remapping

def remap(ref_map, temp_map, out_shape = None, verbose = False):
    """
    reproject hmi map onto hrt with hrt pixel size and observer coordinates
    from SO/PHI-HRT pipeline
    Parameters
    ----------
    ref_map : sunpy.map.GenericMap
        reference map
    temp_map : sunpy.map.GenericMap
        map to be remapped
    out_shape : tuple, None
        shape of output map, if None the out_shape is the shape of ref_map (Default: None)

    Returns
    -------
    temp_remap : sunpy.map.GenericMap
        remapped hmi map
    """
    import sunpy.map
    from reproject import reproject_adaptive
    from sunpy.coordinates import propagate_with_solar_surface
    from astropy.wcs import WCS

    if out_shape is None:
        out_shape = ref_map.data.shape
    
    # define new header for hmi map using hrt observer coordinates
    with propagate_with_solar_surface():
        out_header = sunpy.map.make_fitswcs_header(
            out_shape,
            ref_map.reference_coordinate.replicate(rsun=temp_map.reference_coordinate.rsun),
            scale=u.Quantity(ref_map.scale),
            instrument=temp_map.instrument,
            observatory=temp_map.observatory,
            wavelength=temp_map.wavelength
        )

        out_header['dsun_obs'] = ref_map.coordinate_frame.observer.radius.to(u.m).value
        out_header['hglt_obs'] = ref_map.coordinate_frame.observer.lat.value
        out_header['hgln_obs'] = ref_map.coordinate_frame.observer.lon.value
        out_header['detector'] = temp_map.detector
        out_header['crpix1'] = ref_map.fits_header['CRPIX1']
        out_header['crpix2'] = ref_map.fits_header['CRPIX2']
        out_header['crval1'] = ref_map.fits_header['CRVAL1']
        out_header['crval2'] = ref_map.fits_header['CRVAL2']
        
        if 'CROTA' in ref_map.fits_header:
            key = 'CROTA'
        else:
            key = 'CROTA2'
        out_header[key] = ref_map.fits_header[key]
        
        if 'PC1_1' in ref_map.fits_header:
            out_header['PC1_1'] = ref_map.fits_header['PC1_1']
            out_header['PC1_2'] = ref_map.fits_header['PC1_2']
            out_header['PC2_1'] = ref_map.fits_header['PC2_1']
            out_header['PC2_2'] = ref_map.fits_header['PC2_2']
        else:
            out_header['PC1_1'] = np.cos(ref_map.fits_header[key]*np.pi/180)
            out_header['PC1_2'] = -np.sin(ref_map.fits_header[key]*np.pi/180)
            out_header['PC2_1'] = np.sin(ref_map.fits_header[key]*np.pi/180)
            out_header['PC2_2'] = np.cos(ref_map.fits_header[key]*np.pi/180)

        out_wcs = WCS(out_header)
        
        # reprojection
        temp_origin = temp_map
        output, _ = reproject_adaptive(temp_origin, out_wcs, out_shape,kernel='Hann',boundary_mode='ignore')
        temp_remap = sunpy.map.Map((output, out_header))
    temp_remap.plot_settings = temp_origin.plot_settings

    
    return temp_remap

###############################################
# Function to compute the shift with sub-pixel accuracy
 
def circular_mask(h, w, center, radius):
    """create a circular mask

    Parameters
    ----------
    h : int
        height of the mask
    w : int
        width of the mask
    center : [x,y]
        center of the mask
    radius : float
        radius of the mask

    Returns
    -------
    mask: 2D array
        mask with 1 inside the circle and 0 outside
    """
    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = dist_from_center <= radius
    return mask

def image_register(ref,im,subpixel=True,deriv=False,d=50,verbose=False):
    """
    Computes shift between two maps with sub-pixel accuracy by fitting a gaussian on the cross-correlation map
    from SO/PHI-HRT pipeline

    Parameters
    ----------
    ref: numpy.array
        reference image
    im: numpy.array
        image to be registered
    subpixel: bool
        if True, the shift is computed with sub-pixel accuracy. Default: True
    deriv: bool
        if True, the derivative of the images are used to compute the shift. Default: False
    d: int
        radius of the mask used to fit the gaussian when subpixel is True. It is updated if the fit fails. Default: 50
    verbose: bool
        if True, print information about fitting process. Default: False
    
    Returns
    -------
    r: numpy.array
        cross-correlation map
    shifts: list
        [y,x] shifts etween the images in pixel units
    """
    try:
        import pyfftw.interfaces.numpy_fft as fft
    except:
        import numpy.fft as fft
        
    def _image_derivative(d):
        import numpy as np
        from scipy.signal import convolve
        kx = np.asarray([[1,0,-1], [1,0,-1], [1,0,-1]])
        ky = np.asarray([[1,1,1], [0,0,0], [-1,-1,-1]])
        kx=kx/3.
        ky=ky/3.
        SX = convolve(d, kx,mode='same')
        SY = convolve(d, ky,mode='same')
        A=SX**2+SY**2
        return A

    def _g2d(X, offset, amplitude, sigma_x, sigma_y, xo, yo, theta):
        import numpy as np
        (x, y) = X
        xo = float(xo)
        yo = float(yo)
        a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
        b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
        c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
        g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo)
                                + c*((y-yo)**2)))
        return g.ravel()

    def _gauss2dfit(a,mask):
        import numpy as np
        from scipy.optimize import curve_fit
        sz = np.shape(a)
        X,Y = np.meshgrid(np.arange(sz[1])-sz[1]//2,np.arange(sz[0])-sz[0]//2)
        c = np.unravel_index(a.argmax(),sz)
        Xf = X[mask>0]; Yf = Y[mask>0]; af = a[mask>0]

        stdx = .5 
        stdy = .5 
        initial_guess = [np.median(a), np.max(a), stdx, stdy, c[1] - sz[1]//2, c[0] - sz[0]//2, 0]
        bounds = ([-1,-1,0,0,initial_guess[4]-1,initial_guess[5]-1,-180],
                  [1,1,initial_guess[2]*4,initial_guess[2]*4,initial_guess[4]+1,initial_guess[5]+1,180])

        popt, pcov = curve_fit(_g2d, (Xf, Yf), af.ravel(), p0=initial_guess,bounds=bounds)
        return np.reshape(_g2d((X,Y), *popt), sz), popt
    
    def _one_power(array):
        return array/np.sqrt((np.abs(array)**2).mean())

    if deriv:
        ref = _image_derivative(ref - np.mean(ref))
        im = _image_derivative(im - np.mean(im))
        
    shifts=np.zeros(2)
    FT1=fft.fftn(ref - np.mean(ref))
    FT2=fft.fftn(im - np.mean(im))
    ss=np.shape(ref)
    r=np.real(fft.ifftn(_one_power(FT1) * _one_power(FT2.conj())))
    r = fft.fftshift(r)
    rmax=np.max(r)
    ppp = np.unravel_index(np.argmax(r),ss)
    shifts = [ppp[0]-ss[0]//2,ppp[1]-ss[1]//2]
    if subpixel:
        dd = [d,75,100,30]
        for d1 in dd:
            try:
                if d1>0:
                    mask = circular_mask(ss[0],ss[1],[ppp[1],ppp[0]],d1)
                else:
                    mask = np.ones(ss,dtype=bool); d = ss[0]//2
                g, A = _gauss2dfit(r,mask)
                break
            except RuntimeError as e:
                if verbose: print(f"Issue with gaussian fitting using mask with radius {d1}\nTrying new value...")
                if d1 == dd[-1]:
                    raise RuntimeError(e)
        shifts[0] = A[5]
        shifts[1] = A[4]
        del g
    del FT1, FT2
    return r, shifts

###############################################
# Function to translate the header CRPIX or CRVAL by a certain amount of pixels

def translate_header(h,tvec,mode='crpix'):
    """calculate new header when image is translated by a fixed vector

    Parameters
    ----------
    h : astropy.io.fits.header.Header
        header of image to be translated
    tvec : list
        vector to translate image by (in pixels) [y,x]
    mode : str
        if 'crpix' (Default) the shift will be applied to CRPIX*, if 'crval' the shift will be applied to CRVAL*

    Returns
    -------
    h: astropy.io.fits.header.Header
        new header
    """
    if mode == 'crval':
        tr = np.asarray([[1,0,-tvec[1]],[0,1,-tvec[0]],[0,0,1]])
        if 'CROTA' in h:
            angle = h['CROTA'] # positive angle = clock-wise rotation of the reference system axes 
        else:
            angle = h['CROTA2']
        rad = angle * np.pi/180
        vec = np.asarray([tvec[1],tvec[0],1])
        rot = np.asarray([[np.cos(rad),-np.sin(rad),0],[np.sin(rad),np.cos(rad),0],[0,0,1]])
        shift = rot @ vec
        shift[0] *= h['CDELT1']
        shift[1] *= h['CDELT2']
        h['CRVAL1'] = round(h['CRVAL1']-shift[0],4)
        h['CRVAL2'] = round(h['CRVAL2']-shift[1],4)
    elif mode == 'crpix':
        tr = np.asarray([[1,0,tvec[1]],[0,1,tvec[0]],[0,0,1]])
        coords = np.asarray([h['CRPIX1'],h['CRPIX2'],1])
        new_coords = tr @ coords
        h['CRPIX1'] = round(new_coords[0],4)
        h['CRPIX2'] = round(new_coords[1],4)
    else:
        print('mode not valid\nreturn old header')
    return h

###############################################
# Function to shift an image by applying a phase shift in Fourier space

def fft_shift(img,shift):
    """Shift an image in the Fourier domain and return the shifted image (non fourier domain)

    Parameters
    ----------
    img : 2D-image
        2D-image to be shifted
    shift : list
        [dy,dx] shift in pixel

    Returns
    -------
    img_shf : 2D-image
        shifted image
    """
    try:
        import pyfftw.interfaces.numpy_fft as fft
    except:
        import numpy.fft as fft
    sz = img.shape
    ky = fft.ifftshift(np.linspace(-np.fix(sz[0]/2),np.ceil(sz[0]/2)-1,sz[0]))
    kx = fft.ifftshift(np.linspace(-np.fix(sz[1]/2),np.ceil(sz[1]/2)-1,sz[1]))

    img_fft = fft.fft2(img)
    shf = np.exp(-2j*np.pi*(ky[:,np.newaxis]*shift[0]/sz[0]+kx[np.newaxis]*shift[1]/sz[1]))
    
    img_fft *= shf
    img_shf = fft.ifft2(img_fft).real
    
    return img_shf
    
###############################################
# WCS correction
# --------------------------------------------
# This is an iterative procedure in which the remapped HMI and HRT are cross-correlated to sub-pixel levels.
# Once the shift is found, the HRT header is updated by the shift and HMI is reprojected on the new WCS set.
# This new map is then correlated again to find the new shift. This whole procedure runs until the shift is smaller than a given threshold.

###############################################
# First step, we rotate the SDO/HMI dataset to the same CROTA as SO/PHI-HRT

hmi_map_rot = hmi_map.rotate((hmi_map.fits_header['CROTA2']-hrt_map.fits_header['CROTA'])*u.deg,method='opencv')

###############################################
# Then, we select a region of the SO/PHI-HRT FoV, big enough to avoid the feature to be out of the FoV. You can also choose the whole FoV.

fov = (slice(512,1536),slice(512,1536)) # (y,x)
ht = hrt_map.fits_header
hrt_data = hrt_map.data.copy()

hrt_map = hrt_map.submap(np.asarray([fov[1].start, fov[0].start])*u.pix,
                                    top_right=np.asarray([fov[1].stop-1, fov[0].stop-1])*u.pix)

###############################################
# Iterative procedure

from sunpy.coordinates import Helioprojective
from sunpy.coordinates import propagate_with_solar_surface

t0 = ht['DATE-AVG']
t_obs = datetime.datetime.fromisoformat(ht['DATE-AVG'])

h_hrt = hrt_map.fits_header.copy()
ht['DATE-BEG'] = t0; ht['DATE-OBS'] = t0

shift = [1,1]
i = 0
match = True
max_iterations = 10
THRESHOLD_REMAP = 5E-2
THRESHOLD_SHIFT = 1E-2


h_rot = hmi_map_rot.fits_header
h_rot['CROTA2'] = hrt_map.fits_header['CROTA']

hmi_map_rot = sunpy.map.Map((np.nan_to_num(hmi_map_rot.data,nan=0,posinf=0,neginf=0),h_rot))

while np.any(np.abs(shift)>THRESHOLD_REMAP):
    hrt_map = sunpy.map.Map((hrt_map.data.copy(),h_hrt))
    
    with propagate_with_solar_surface():
        with SphericalScreen(hmi_map_rot.observer_coordinate, only_off_disk=True):
            
            bl = hrt_map.bottom_left_coord
            tr = hrt_map.top_right_coord

            top_right = hmi_map_rot.world_to_pixel(tr)
            bottom_left = hmi_map_rot.world_to_pixel(bl)

            tr_hmi_map = np.array([top_right.x.value,top_right.y.value])
            bl_hmi_map = np.array([bottom_left.x.value,bottom_left.y.value])

            hmi_submap = hmi_map_rot.submap(bl_hmi_map*u.pix,top_right=tr_hmi_map*u.pix)
            # hmi_submap.fits_header['CROTA2'] = hmi_map_rot.fits_header['CROTA2']

            hpc_coords = sunpy.map.all_coordinates_from_map(hmi_submap)
            mask = ~sunpy.map.coordinate_is_on_solar_disk(hpc_coords)

            tt = hmi_submap.data; tt[mask==1] = 0
            hmi_submap = sunpy.map.Map((tt,hmi_submap.fits_header)); del tt

            hrt_remap = remap(hmi_submap, hrt_map)
                
            ref = np.nan_to_num(hmi_submap.data.copy(),True,0,0,0)
            temp = np.nan_to_num(hrt_remap.data.copy(),True,0,0,0)

    s = [1,1]
    shift = [0,0]
    it = 0

    deriv = False

    while np.any(np.abs(s)>THRESHOLD_SHIFT) and it<10:
        if it == 0:
            _,s = image_register(ref,temp,False,deriv)
            if np.any(np.abs(s)==0):
                _,s = image_register(ref,temp,True,deriv)
        else:
            _,s = image_register(ref,temp,True,deriv)
        
        shift = [shift[0]+s[0],shift[1]+s[1]]
        temp = np.nan_to_num(fft_shift(np.nan_to_num(hrt_remap.data.copy(),True,0,0,0), shift),True,0,0,0)
        it += 1

    ht = translate_header(ht.copy(),-np.asarray(shift)*hmi_submap.fits_header['CDELT1']/hrt_map.fits_header['CDELT1'] * \
                                    hmi_map_rot.fits_header['DSUN_OBS']/hrt_map.fits_header['DSUN_OBS'],
                                    mode='crval')
    h_hrt = translate_header(h_hrt.copy(),-np.asarray(shift)*hmi_submap.fits_header['CDELT1']/hrt_map.fits_header['CDELT1'] * \
                                    hmi_map_rot.fits_header['DSUN_OBS']/hrt_map.fits_header['DSUN_OBS'],
                                    mode='crval')
    
    print(it,'iterations shift (x,y):',round(shift[1],2),round(shift[0],2))

    i+=1
    if i == max_iterations:
        print('Maximum iterations reached:',i)
        match = False
        break

hrt_map = sunpy.map.Map((hrt_data,ht))

###############################################
# Check the result
# --------------------------------------------
from matplotlib.widgets import Slider
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Create the figure and axes
fig, ax = plt.subplots(figsize=(8,8))
plt.subplots_adjust(top=0.85)

# Display the first image
img1 = ax.imshow(hrt_map.data, cmap='hmimag', clim=(-1500,1500), origin='lower')
hmi_remap = hmi_map_rot.reproject_to(hrt_map.wcs)
# Display the second image with initial alpha
img2 = ax.imshow(hmi_remap.data, cmap='hmimag', clim=(-1500,1500), origin='lower', alpha=0.5)

# Create slider axis and slider
slider_ax = inset_axes(ax, width="60%", height="5%", loc='upper center', borderpad=0.6)
slider = Slider(slider_ax, '', 0.0, 1.0, valinit=0.5)
label_text = slider_ax.text(0.5, 1.6, 'Transparency [HRT --> HMI]', transform=slider_ax.transAxes,
                            ha='center', va='bottom', fontsize=10)

# Update function for the slider
def update(val):
    alpha = slider.val
    img2.set_alpha(alpha)
    fig.canvas.draw_idle()

slider.on_changed(update)

plt.show()


###############################################
# Print Values
# --------------------------------------------

keys = ['PHIDATID', 'CROTA', 'PC1_1', 'PC1_2', 'PC2_1', 'PC2_2', 'CRPIX1', 'CRPIX2', 'CRVAL1', 'CRVAL2']
ends = [',']*(len(keys)-1)+['\n']
cols = ''
for e,k in zip(ends,keys):
    cols += k+e

cols_old = ''
for e,k in zip(ends,keys):
    cols_old += k+e

for e,k in zip(ends,keys):
        v = hrt_map.fits_header[k]
        if isinstance(v,str):
            cols += v+e
        else:
            cols += '{:{width}.{prec}f}'.format(v,width=5,prec=3)+e
        
        v = original_header[k]
        if isinstance(v,str):
            cols_old += v+e
        else:
            cols_old += '{:{width}.{prec}f}'.format(v,width=5,prec=3)+e

print('New Header')
print(cols)
print('')
print('Old Header')
print(cols_old)