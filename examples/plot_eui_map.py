"""
==================================
Finding and Plotting EUI Data
==================================

This example demonstrates how to:

* Search for Solar Orbiter EUI (Extreme Ultraviolet Imager) observations
* Download a sample image
* Plot it using SunPy

"""

import matplotlib.pyplot as plt
import sunpy.map
from sunpy.net import Fido, attrs as a
import sunpy_soar

###############################################################################
# Registering the Solar Orbiter Archive (SOAR)
# --------------------------------------------
#
# By importing `sunpy_soar`, the Solar Orbiter Archive (SOAR) is **automatically**
# registered as a data provider in `Fido`. This allows us to directly query and 
# download Solar Orbiter data just like any other SunPy Fido search.
#
# Now, we will use `Fido` to search for available EUI Full Sun Imager (FSI) 174 angstrom data, namely
# the dataproduct called **EUI-FSI174-IMAGE**, within a specific time range.
#
# The search will return metadata about available files.

search_results = Fido.search(a.Time("2022-03-01", "2022-03-01 01:00"),
                             a.soar.Product("EUI-FSI174-IMAGE"),
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
# Visualizing the EUI Image with `sunpy.map`
# -------------------------------------------
#
# The downloaded file is in **FITS format**, which is commonly used for astronomical
# imaging data. We use `sunpy.map.Map` to load and display it.

eui_map = sunpy.map.Map(downloaded_files[0])

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(projection=eui_map)
eui_map.plot(axes=ax)
plt.colorbar()
plt.title("Solar Orbiter EUI Image")

plt.show()
