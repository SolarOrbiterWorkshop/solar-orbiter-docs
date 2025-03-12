"""
==================================
Finding and Plotting Solar Orbiter EUI Data
==================================

This example demonstrates how to:
- Search for Solar Orbiter EUI (Extreme Ultraviolet Imager) data
- Download a sample image
- Plot it using SunPy
"""

import matplotlib.pyplot as plt
import sunpy.map
from sunpy.net import Fido, attrs as a
import sunpy_soar

###############################################################################
# First, we will use `Fido <sunpy.net.fido_factory.UnifiedDownloaderFactory>`__
# to search for available Solar Orbiter **EUI-FSI174-IMAGE** data in a specific
# time range.  
#
# The search will return metadata about available files.

search_results = Fido.search(a.Time("2022-03-01", "2022-03-01 01:00"),
                             a.soar.Product("EUI-FSI174-IMAGE"),
                             a.Level(2))

print(search_results)

###############################################################################
# Once we have the search results, we can fetch the **first available file**.

downloaded_files = Fido.fetch(search_results[0, 0])

###############################################################################
# Now that we have the file, we can create a `~sunpy.map.Map` to visualize it.

eui_map = sunpy.map.Map(downloaded_files[0])

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(projection=eui_map)
eui_map.plot(axes=ax)
plt.colorbar()
plt.title("Solar Orbiter EUI Image")

plt.show()
