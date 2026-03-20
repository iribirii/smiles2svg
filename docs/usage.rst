Usage
=====

Command line
------------

After installation, the ``smiles2svg`` script is on your ``PATH``:

.. code-block:: bash

   smiles2svg -s 'Cn1cnc2n(C)c(=O)n(C)c(=O)c12'

You can also run the package as a module:

.. code-block:: bash

   python -m smiles2svg -s 'c1ccccc1'

See all options:

.. code-block:: bash

   smiles2svg -h

Batch input: one SMILES per line in a text file:

.. code-block:: bash

   smiles2svg -f examples/smiles.txt

PNG export (requires Cairo):

.. code-block:: bash

   smiles2svg -s 'CCO' --png --png-width 900

Canvas (xyzrender-style)
--------------------------

The structure is laid out in **pixel space**: the longer in-plane axis of the molecule is scaled to fit ``canvas_size − 2×padding`` (default 800 px and 20 px padding). Let ``scale`` = (that inner size in px) / (molecular bounding-box span in Å). Bond stroke width and label ``font-size`` are ``constant_Å × scale`` (defaults match legacy 0.1 Å / 0.8 Å stroke and font in coordinate space), so **proportions stay the same** when you change canvas size or when the molecule is larger or smaller in Å. The `xyzrender`_ project uses a related pixel pipeline; stroke/font sizing there is configured differently but the idea of a uniform ``scale`` is the same.

.. _xyzrender: https://github.com/aligfellow/xyzrender

.. code-block:: bash

   smiles2svg -s 'CCO' --canvas-size 800 --padding 20

Python API
----------

.. code-block:: python

   from smiles2svg import render_smiles_to_svg

   svg = render_smiles_to_svg(
       "CCO",
       style="plain",
       color="default",
       canvas_size=800,
       padding=20.0,
       # output_path="molecule.svg",  # optional: write to disk
   )
   print(svg[:200])

Drawing styles
--------------

``--style`` / ``style`` accepts:

* ``plain`` — colored disks, carbon without labels
* ``names_hetero`` — element symbols on heteroatoms
* ``names_all`` — symbols on every atom
* ``stroke`` — colored disks with outline

See the project README for example images.
