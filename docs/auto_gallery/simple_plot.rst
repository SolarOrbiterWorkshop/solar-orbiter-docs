
.. DO NOT EDIT.
.. THIS FILE WAS AUTOMATICALLY GENERATED BY SPHINX-GALLERY.
.. TO MAKE CHANGES, EDIT THE SOURCE PYTHON FILE:
.. "auto_gallery/simple_plot.py"
.. LINE NUMBERS ARE GIVEN BELOW.

.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        :ref:`Go to the end <sphx_glr_download_auto_gallery_simple_plot.py>`
        to download the full example code.

.. rst-class:: sphx-glr-example-title

.. _sphx_glr_auto_gallery_simple_plot.py:


Simple Plot Example
===================

This example demonstrates a basic plot of solar X-ray flux.

.. GENERATED FROM PYTHON SOURCE LINES 8-21

.. code-block:: Python


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


.. _sphx_glr_download_auto_gallery_simple_plot.py:

.. only:: html

  .. container:: sphx-glr-footer sphx-glr-footer-example

    .. container:: sphx-glr-download sphx-glr-download-jupyter

      :download:`Download Jupyter notebook: simple_plot.ipynb <simple_plot.ipynb>`

    .. container:: sphx-glr-download sphx-glr-download-python

      :download:`Download Python source code: simple_plot.py <simple_plot.py>`

    .. container:: sphx-glr-download sphx-glr-download-zip

      :download:`Download zipped: simple_plot.zip <simple_plot.zip>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
