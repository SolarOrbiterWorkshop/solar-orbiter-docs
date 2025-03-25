"""
=====================================
Plotting the Solar Orbiter Trajectory
=====================================

This example demonstrates how to to retrieve Solar Orbiter's position using SunPy and JPL HORIZONS, through :func:`sunpy.coordinates.get_horizons_coord`.
We will get the position of Solar Orbiter as a function of time and plot this in a polar plot of lon and radius in heliographic stonyhurst frame.
"""

import matplotlib.pyplot as plt
from astropy.time import Time
import astropy.units as u
from sunpy.coordinates import get_horizons_coord

###############################################################################
# Get Solar Orbiter’s Ephemeris
# -----------------------------
#
# We use SunPy’s :func:`sunpy.coordinates.get_horizons_coord` to query Solar Orbiter's position
# from JPL HORIZONS. This returns coordinates in a SkyCoord object, in the Heliographic Stonyhurst frame.

times = Time("2024-01-01") + (u.Quantity(range(0, 180, 1)) * u.day)
solo_coords = get_horizons_coord("Solar Orbiter", times)

###############################################################################
# Plotting in the Heliographic Stonyhurst (HGS) Frame
# ---------------------------------------------------
#
# We plot Solar Orbiter's radial distance from the Sun against its
# heliographic longitude using a polar plot. This gives a clear view
# of the spacecraft's orbit around the Sun.

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(projection="polar")

ax.plot(solo_coords.lon.to(u.rad), 
        solo_coords.radius.to(u.AU),
        label='Solar Orbiter', marker='.', ls='-')

# Mark the Sun at the origin
ax.scatter(0, 0, marker='o', color='orange', label="Sun")

ax.set_theta_zero_location('S')  
ax.set_theta_direction(-1)       
ax.set_title('Solar Orbiter Orbit in Heliographic Stonyhurst (HGS) from {:s} to {:s}'.format(times[0].strftime("%Y-%m-%d"),times[0].strftime("%Y-%m-%d")), 
	         pad=20)
ax.legend()

plt.tight_layout()
plt.show()
