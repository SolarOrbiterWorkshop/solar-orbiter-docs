==================================
Contributing to the Documentation
==================================

Thank you for your interest in improving the **Solar Orbiter Data Documentation**!

There are several ways to contribute:

- Update documentation: Fix typos, improve explanations, or add useful links.
- Add gallery examples: Contribute Python scripts demonstrating Solar Orbiter data analysis.
- Request new examples: Suggest examples you’d like to see.
- Report issues: Found a problem or missing information? Open an issue on GitHub.


Requesting an Example or Reporting an Issue
===========================================

If you have an idea for a new example or a problem with the documentation,  
**open an issue on GitHub**: `GitHub Issues <https://github.com/SolarOrbiterWorkshop/solar-orbiter-docs/issues>`__

When submitting an issue, please include:
 * A short description of the example you’d like to see.
 * If possible, links to related data or papers.



Getting Started with Contributions
==================================

This documentation is hosted on GitHub.  
To make changes, follow these steps:

1. Fork the repository on GitHub:  
   `Solar Orbiter Docs Repository <https://github.com/SolarOrbiterWorkshop/solar-orbiter-docs>`__
2. Clone your fork:

   .. code-block:: bash

      git clone https://github.com/YOUR_USERNAME/solar-orbiter-docs.git
      cd solar-orbiter-docs

3. Create a new branch for your changes:

   .. code-block:: bash

      git checkout -b update-docs

4. Make edits to the `docs/` folder (e.g., update `.rst` files, add gallery examples).
5. Test your changes locally:

   .. code-block:: bash

      cd docs
      make clean
      make html

   Open `_build/html/index.html` in a browser to check your updates.
6. Commit and push your changes:

   .. code-block:: bash

      git add .
      git commit -m "Update documentation"
      git push origin update-docs

7. Open a Pull Request on GitHub.

For guidance on using GitHub, see the `Creating a pull request guide <https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request>`__.


Adding an Example to the Gallery
===================================

To contribute an example to the **Example Gallery**, follow these steps:

1. Create a Python script in `examples/` (e.g., `examples/my_example.py`).
2. The script should include:

   - A **docstring** explaining its purpose.

   - **A plot or output** using `matplotlib`, `sunpy`, or other relevant libraries.

3. Example format:

   .. code-block:: python

      """
      Example: Plotting Solar Orbiter EUI Data
      """
      import sunpy.map
      import matplotlib.pyplot as plt

      # Load example EUI data
      my_map = sunpy.map.Map("path/to/eui_data.fits")

      # Plot the map
      my_map.plot()
      plt.show()

4. Add your script to the gallery:

   * Place it inside `examples/`
   * It will automatically appear in the **Example Gallery** after merging.



Adding Useful Links
===================================

If you know of useful Solar Orbiter-related sites, add them to `index.rst` under the **Useful Links** section.

Example:

   .. code-block:: rst

      `Solar Orbiter Science Nuggets <https://www.cosmos.esa.int/web/solar-orbiter/science-nuggets>`__

---


Technical Details
===================================

- This documentation is built using **Sphinx** and **PyData Sphinx Theme**.
- Builds are automatically deployed via **Read the Docs**.
- To install dependencies locally, run:

   .. code-block:: bash

      pip install -r requirements.txt
