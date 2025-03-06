"""
Finding and Plotting Solar Orbiter EUI Data
===========================================

This example demonstrates how to:
- Search for Solar Orbiter EUI (Extreme Ultraviolet Imager) data
- Download a sample image
- Plot it using SunPy

"""

import sunpy.map
import matplotlib.pyplot as plt
from sunpy.net import Fido, attrs as a
import sunpy_soar

# Search for Solar Orbiter EUI images in a given time range
search_results = Fido.search(a.Time("2022-03-01", "2022-03-01 01:00"), 
				             a.soar.Product("EUI-FSI174-IMAGE"), 
				             a.Level(2))

# Print search results
print(search_results)

# Fetch the first available file
downloaded_files = Fido.fetch(search_results[0, 0])

# Load the downloaded EUI image
eui_map = sunpy.map.Map(downloaded_files[0])

# Plot the EUI image
plt.figure(figsize=(8, 6))
eui_map.plot()
plt.colorbar()
plt.title("Solar Orbiter EUI Image")
plt.show()
