Getting Started
===============

This page provides an overview of how to set up and use the Solar Orbiter data tools.

Installation
------------

To install required dependencies, run:

.. code-block:: bash

    pip install sunpy sunpy_soar numpy matplotlib

Basic Example
-------------

Here's how to load and plot some Solar Orbiter data:

.. code-block:: python

    import sunpy.map
    from sunpy.net import Fido, attrs as a

    result = Fido.search(a.Time("2024-03-01", "2024-03-02"), a.Instrument.eui)
    print(result)
