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

Install from PyPI
------------------

.. code-block:: bash

   pip install smiles2svg

Install from source
-------------------

.. code-block:: bash

   git clone https://github.com/iribirii/smiles2svg.git
   cd smiles2svg
   pip install .

Using `uv <https://docs.astral.sh/uv/>`_:

.. code-block:: bash

   uv pip install -e ".[dev,docs]"

Documentation build (Furo + Sphinx):

.. code-block:: bash

   pip install -e ".[docs]"
   cd docs && make html

Regenerate documentation figures (optional):

.. code-block:: bash

   python docs/scripts/generate_doc_figures.py
