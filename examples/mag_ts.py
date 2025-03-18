"""
=========================================
Finding and Plotting Solar Orbiter MAG Data
=========================================

This example demonstrates how to search for, download, and visualize 1-minute averaged Solar Orbiter 
MAG (Magnetometer) data from the SOAR archive, load it into a SunPy TimeSeries, and plot it.
"""

import matplotlib.pyplot as plt
from sunpy.net import Fido, attrs as a
import sunpy_soar
from sunpy.timeseries import TimeSeries as ts

###############################################################################
#
# By importing `sunpy_soar`, the Solar Orbiter Archive (SOAR) is **automatically**
# registered as a data provider in `Fido`. This allows us to directly query and 
# download Solar Orbiter data just like any other SunPy Fido search.
#
# Now, we will use `Fido` to search for MAG Level 2 data in Radial-Tangential-Normal (RTN) coordinates.
# We will request 1-minute averaged data over a two-day period.

search_results = Fido.search(
    a.Time("2022-09-04", "2022-09-06 00:00"),
    a.Instrument.mag,
    a.soar.Product('mag-rtn-normal-1-minute'),
    a.Level(2)  
)

print(search_results)

###############################################################################
# The search results contain a list of available files matching our query.
# We use `Fido.fetch` to download the files.

downloaded_files = Fido.fetch(search_results)

###############################################################################
# Loading the MAG Data into a SunPy TimeSeries
# -------------------------------------------
#
# The downloaded file contains time-series magnetic field measurements.
# We use `sunpy.timeseries.TimeSeries` to load and analyse the data.
#
# The **concatenate=True** argument ensures that if multiple files are downloaded,
# they are combined into a single continuous TimeSeries object.

mag_ts = ts(downloaded_files, concatenate=True)

###############################################################################
# Plotting the Magnetic Field Components
# --------------------------------------
# We now plot these components over time.

mag_ts.plot(columns=["B_RTN_0", "B_RTN_1", "B_RTN_2"])

plt.title("Solar Orbiter MAG Data (RTN Coordinates)")
plt.show()
