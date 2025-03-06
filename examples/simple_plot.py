"""
Simple Plot Example
===================

This example demonstrates a basic plot of solar X-ray flux.

"""

import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0, 10, 100)
y = np.sin(x)

plt.figure()
plt.plot(x, y)
plt.xlabel("Time (s)")
plt.ylabel("Flux")
plt.title("Solar X-ray Flux Example")
plt.show()
