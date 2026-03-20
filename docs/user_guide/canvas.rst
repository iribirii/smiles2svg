Pixel canvas and proportions
===============================

Coordinates from RDKit are in **Å**. The renderer maps them to **SVG pixel space**:

* ``scale`` = (``canvas_size`` − 2×``padding``) / max(bounding-box span in Å)
* Bond stroke width and label font size scale as **fixed factors × scale** (defaults match legacy 0.1 Å / 0.8 Å in coordinate space), so **relative** thickness stays consistent when you change ``--canvas-size`` or draw larger/smaller molecules.

Defaults match the **xyzrender**-style idea of an explicit ``viewBox="0 0 W H"`` with ``width``/``height`` in pixels.

See the ``layout`` module in the :doc:`../api` reference for implementation details.
