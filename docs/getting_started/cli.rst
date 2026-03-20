Command-line interface
=======================

Run ``smiles2svg -h`` for the full help text. Summary:

Input
-----

``-s / --smiles`` — one SMILES string.

``-f / --smiles-file`` — path to a text file with **one SMILES per line** (cannot be combined with ``-s``).

``-n / --name`` — output basename (``.svg`` / ``.png``) when a **single** ``-s`` is used. Ignored for batch ``-f``.

Drawing
-------

``--style`` — ``plain`` | ``names_hetero`` | ``names_all`` | ``stroke``

``--color`` — hex color for all atoms, or ``default`` for per-element colors.

``--bond-color`` — bond stroke (default ``#717171``).

``--font`` — CSS ``font-family`` stack for labels (names\_\* styles).

``--add-hydrogens`` — draw implicit H on typical H-bond donors (O, N, P).

Output
------

``--png`` — also write PNG (needs Cairo).

``--png-width`` — PNG width in pixels (default 900).

``--canvas-size`` — SVG canvas target (default 800 px; longer axis of the molecule is fitted).

``--padding`` — inner padding in px (default 20).

Example
-------

.. code-block:: bash

   smiles2svg -s 'CCO' --style names_hetero --png -n ethanol
