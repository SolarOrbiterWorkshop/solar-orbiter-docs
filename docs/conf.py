import os
import sys
from datetime import datetime

# Add project root to path
sys.path.insert(0, os.path.abspath('..'))

# Project information
project = "Solar Orbiter Tools"
author = "Laura Hayes"
copyright = f"{datetime.now().year}, {author}"
release = "0.1"

# -- General configuration ---------------------------------------------------

source_suffix = ".rst"

# The master toctree document.
master_doc = "index"

# Sphinx extensions (including Sphinx-Gallery)
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_gallery.gen_gallery",  # Example gallery
    "sphinx_design",  # For the grid-style landing page
    "sphinx.ext.intersphinx",
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "astropy": ("https://docs.astropy.org/en/stable/", None),
    "sunpy": ("https://docs.sunpy.org/en/stable/", None),
}

# Paths
templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output -------------------------------------------------
# Static files
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_logo = "_static/SolarOrbiter.png"
html_theme = "pydata_sphinx_theme"


# Sidebar and Navigation
# html_sidebars = {
#     "index": ["search-button-field"],
#     "**": ["search-button-field", "sidebar-nav-bs"]
# }

# Custom Theme Options
html_theme_options = {
    "logo": {
        "text": "Solar Orbiter Docs",
    },
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/SolarOrbiterWorkshop/solar-orbiter-docs",
            "icon": "fa-brands fa-github",
        },
        {
            "name": "ESA Solar Orbiter",
            "url": "https://www.esa.int/Science_Exploration/Space_Science/Solar_Orbiter",
            "icon": "fa-solid fa-sun",
        },
        {
            "name": "SunPy",
            "url": "https://docs.sunpy.org/",
            "icon": "fa-brands fa-python",
        },
    ],
    "navbar_start": ["navbar-logo"],
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "secondary_sidebar_items": ["page-toc"],
}



# -- Sphinx-Gallery Configuration (Example Gallery) --------------------------
os.environ["JSOC_EMAIL"] = "jsoc@sunpy.org"
sphinx_gallery_conf = {
    'examples_dirs': '../examples',   # Path to example scripts
    'gallery_dirs': 'auto_gallery',   # Output path for gallery
    'filename_pattern': r'.*\.py',    # Only include .py files
    'download_all_examples': False,   # Don't zip all examples for download
    'remove_config_comments': True,   # Clean up script headers
}

# Source file suffix
source_suffix = ".rst"
master_doc = "index"

# Custom setup
def setup(app):
    app.add_css_file("custom.css")  # Add custom CSS for styling
