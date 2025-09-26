"""
===============================================
Finding and Plotting Solar Orbiter EPD/EPT Data
===============================================

This example demonstrates how to search for, download, and visualize Level-3
1-minute averaged data of the Electron Proton Telescope (EPT) of Solar Orbiter's
Energetic Particle Detector (EPD) from the SOAR archive, load it into a SunPy
TimeSeries, and plot it.
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
# Now, we will use `Fido` to search for EPT Level-3 data.
# We will request 1-minute averaged data over a two-day period.

search_results = Fido.search(
    a.Time("2022-09-04", "2022-09-06 00:00"),
    a.Instrument.epd,
    a.soar.Product('epd-ept-1min'),
    a.Level(3)  
)

print(search_results)

###############################################################################
# The search results contain a list of available files matching our query.
# We use `Fido.fetch` to download the files.

downloaded_files = Fido.fetch(search_results)

###############################################################################
# Loading the EPD/EPT Data into a SunPy TimeSeries
# ------------------------------------------------
#
# The downloaded file contains time-series of energetic charged particle measurements.
# We use `sunpy.timeseries.TimeSeries` to load and analyse the data.
#
# The **concatenate=True** argument ensures that if multiple files are downloaded,
# they are combined into a single continuous TimeSeries object.

ept_ts = ts(downloaded_files, concatenate=True)

###############################################################################
# Plotting Electron Data for Different Viewing Directions
# -------------------------------------------------------
# 
# The EPT instrument measures electron and ion (mostly protons) fluxes in four different viewing directions. 
# They are indicated in the column names as ``S`` (sunward along the nominal Parker spiral),
# ``A`` (anti-sunward), ``N`` (northward), and ``D`` (southward; D stands for "down").
# In addition, the ending of the column names indicates the energy channel (which are **not** the same for electrons and ions).
# Note that the electron data is in the Level-3 data product corrected for contamination.
# The resulting column names are for example ``Electron_Corrected_Flux_S_1`` for the sunward viewing direction and energy channel 1, or ``Ion_Flux_A_3`` for the anti-sunward viewing direction and energy channel 3.
#
# We now plot the different viewing directions of one electron energy channel over time.

channel = 1  # Change this to plot different energy channels 
ept_ts.plot(columns=[f"Electron_Corrected_Flux_{viewing}_{channel}" for viewing in ['S', 'A', 'N', 'D']])

plt.yscale("log")
plt.title("Solar Orbiter EPD/EPT")
plt.show()
