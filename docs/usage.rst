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

Python API
----------

.. code-block:: python

   from smiles2svg import render_smiles_to_svg

   svg = render_smiles_to_svg(
       "CCO",
       style="plain",
       color="default",
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
