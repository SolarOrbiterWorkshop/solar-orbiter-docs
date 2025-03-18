Getting Started with Solar Orbiter Data
========================================

This page provides an overview of how to set up your Python environment and start working with Solar Orbiter data.

This documentation and the Example Gallery build upon existing tools in the SunPy ecosystem wherever possible, leveraging community-developed libraries for solar physics research.


Installing the Required Packages
--------------------------------

Before working with Solar Orbiter data, you need Python and a few key libraries.  
For a detailed guide, including how to install Python, see the `SunPy Installation Guide <https://docs.sunpy.org/en/stable/tutorial/installation.html>`_.

For a quick setup, you can install the core dependencies using pip:

.. code-block:: bash

    pip install sunpy[all] sunpy_soar 

Or, if you're using conda, install them via:

.. code-block:: bash

    conda install -c conda-forge sunpy sunpy-soar 

These packages provide:
- sunpy: Core tools for solar data analysis.
- sunpy_soar: Interface for accessing Solar Orbiter data from the SOAR archive.

Verifying Your Installation
---------------------------

To ensure everything is set up correctly, you can run:

.. code-block:: python

    import sunpy
    import sunpy_soar

    print("SunPy version:", sunpy.__version__)

If this runs without errors, you're ready to start working with Solar Orbiter data!

---

First Example: Searching for Solar Orbiter Data
-----------------------------------------------

Now, let's search for **Solar Orbiter EUI data** from the SOAR archive:

.. code-block:: python

    from sunpy.net import Fido, attrs as a
    import sunpy_soar

    # Search for EUI data within a specific time range
    result = Fido.search(a.Time("2021-02-01 00:00", "2021-02-01 04:00"), a.Instrument.eui, a.Level(2))

    # Display the available data products
    print(result)

This query searches for EUI (Extreme Ultraviolet Imager) data from Solar Orbiter.

Next Steps
----------

Once your setup is complete, check out the following resources:

* :doc:`data_overview/analysis_tools`  
  Learn about additional community-developed Python tools for working with Solar Orbiter data.

* :doc:`auto_gallery/index`  
  Browse practical examples of how to query, download, and analyze Solar Orbiter observations.

* :doc:`contributing`  
  Want to contribute examples or improve the documentation? Find out how you can help!


