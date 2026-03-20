Installation
============

Requirements
------------

* Python 3.9+
* `RDKit <https://www.rdkit.org/>`_ — cheminformatics and 2D coordinates
* `svgwrite <https://svgwrite.readthedocs.io/>`_ — SVG generation
* `NumPy <https://numpy.org/>`_ — geometry helpers
* `CairoSVG <https://cairosvg.org/>`_ — optional PNG export (requires Cairo system libraries)

On Linux, install Cairo development headers (e.g. ``libcairo2-dev`` on Debian/Ubuntu) if you need PNG output.

Install from PyPI (when published)
-----------------------------------

.. code-block:: bash

   pip install smiles2svg

Install from source
-------------------

.. code-block:: bash

   git clone https://github.com/iribirii/smiles2svg.git
   cd smiles2svg
   pip install .

Editable install for development:

.. code-block:: bash

   pip install -e ".[dev,docs]"
