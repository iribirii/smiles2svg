Python API
==========

.. code-block:: python

   from smiles2svg import render_smiles_to_svg

   svg = render_smiles_to_svg(
       "c1ccccc1",
       style="names_hetero",
       color="default",
       bond_color="#717171",
       font="Calibri, Helvetica, Arial, sans-serif",
       add_hydrogens=False,
       canvas_size=800,
       padding=20.0,
       output_path="benzene.svg",  # optional
   )

Invalid SMILES raises ``ValueError``.

Full signatures are documented under :doc:`../api`.
