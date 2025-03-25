===================================
Solar Orbiter Data Overview
===================================

Solar Orbiter is a **European Space Agency (ESA) mission** dedicated to studying the **Sun and the heliosphere**.  
It carries **ten scientific instruments** to observe the Sun's surface, corona, and the surrounding space environment.  
This page provides an **overview of Solar Orbiter data products** and how to access them.


Solar Orbiter Instruments
=========================
Solar Orbiter carries **ten scientific instruments**, classified as **remote-sensing** or **in-situ** instruments.

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Instrument
     - Description
   * - **EUI** (Extreme Ultraviolet Imager)
     - Captures high-resolution images of the solar corona in extreme ultraviolet.
   * - **Metis** (Visible Light and UV Coronagraph)
     - Observes the solar corona in visible and ultraviolet light.
   * - **PHI** (Polarimetric and Helioseismic Imager)
     - Measures the Sun’s magnetic field and helioseismic activity.
   * - **SoloHI** (Heliospheric Imager)
     - Images the solar wind and transient disturbances like coronal mass ejections.
   * - **SPICE** (Spectral Investigation of the Coronal Environment)
     - Analyzes the plasma properties of the solar atmosphere via EUV spectroscopy.
   * - **STIX** (Spectrometer Telescope for Imaging X-rays)
     - Detects X-ray emissions to study solar flares.
   * - **EPD** (Energetic Particle Detector)
     - Measures the properties of energetic particles in the solar wind.
   * - **MAG** (Magnetometer)
     - Measures the magnetic field in the solar wind.
   * - **RPW** (Radio and Plasma Waves)
     - Observes radio emissions and plasma waves in the solar wind.
   * - **SWA** (Solar Wind Analyser)
     - Analyzes the properties of solar wind particles.

More details about each instrument can be found in the `Instrument Documentation <https://www.cosmos.esa.int/web/soar/instrument-documentation>`__.

Solar Orbiter Data Products
===========================
Solar Orbiter data is publicly available via the `Solar Orbiter Archive (SOAR) <https://soar.esac.esa.int/soar/>`_.

Each instrument provides different levels of processed data:

- **Level 0**: Raw instrument data.
- **Level 1**: Calibrated and corrected data.
- **Level 2**: Fully processed, ready-to-use data.


Accessing Solar Orbiter Data
============================
The easiest way to access Solar Orbiter data is through SunPy’s `Fido` search interface.

Example query:

.. code-block:: python

    from sunpy.net import Fido, attrs as a
    import sunpy_soar  # Registers the SOAR data source

    # Search for EUI 174 Å images in a given time range
    results = Fido.search(a.Time("2023-01-01", "2023-01-02"),
                          a.soar.Product("EUI-FSI174-IMAGE"),
                          a.Level(2))

    print(results)

To learn more, see the `Example Gallery <auto_gallery/index.html>`_.


Useful Links
============

- `ESA Solar Orbiter Homepage <https://www.esa.int/Science_Exploration/Space_Science/Solar_Orbiter>`_
- `Solar Orbiter Archive (SOAR) <https://soar.esac.esa.int/soar/>`_
- `Instrument Documentation <https://www.cosmos.esa.int/web/soar/instrument-documentation>`_
- `Solar Orbiter Publications <https://www.cosmos.esa.int/web/solar-orbiter/publications>`_

For contributions, suggestions, or reporting issues, see the :doc:`Contributing Guide </contributing>`.



---

