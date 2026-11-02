# smiles2svg

[![PyPI version](https://img.shields.io/pypi/v/smiles2svg.svg)](https://pypi.org/project/smiles2svg/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Documentation Status](https://readthedocs.org/projects/smiles2svg/badge/?version=latest)](https://smiles2svg.readthedocs.io/en/latest/?badge=latest)

Convert [SMILES](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system) strings into 2D molecular structure diagrams as **SVG**, with optional **PNG** export (via [Cairo](https://www.cairographics.org/)).

**Requires Python 3.10+.**

## Install

```bash
pip install smiles2svg
```

From source (this repository):

```bash
pip install .
```

Runtime dependencies are declared in `pyproject.toml` (`rdkit`, `svgwrite`, `numpy`, `cairosvg`). On Linux, install Cairo development packages (e.g. `libcairo2-dev`) if you need PNG output.

## Command line

```bash
smiles2svg -s 'Cn1cnc2n(C)c(=O)n(C)c(=O)c12'
```

Creates `Cn1cnc2n(C)c(=O)n(C)c(=O)c12.svg` in the current directory.

Rendering uses a **pixel canvas** (similar to [xyzrender](https://github.com/aligfellow/xyzrender)): coordinates are mapped to SVG space with `viewBox="0 0 W H"`, explicit width/height in pixels. Bond stroke and label size are **Å × scale** (defaults 0.10 Å and 0.80 Å for stroke and font, matching legacy proportions), where `scale` is pixels per Å from fitting the molecule to the canvas — so **relative** thickness and text size stay consistent when you change `--canvas-size` or draw a bigger/smaller structure. Defaults: `--canvas-size 800` and `--padding 20`.

```bash
python -m smiles2svg -h
```

### Batch input

A text file with one SMILES per line:

```bash
smiles2svg -f examples/smiles.txt
```

### Custom output basename (single structure only)

```bash
smiles2svg -s 'c1ccccc1' -n benzene
```

Writes `benzene.svg` (and `benzene.png` if `--png` is set).

### PNG export

```bash
smiles2svg -s 'CCO' --png --png-width 900
```

## Python API

```python
from smiles2svg import render_smiles_to_svg

svg = render_smiles_to_svg("c1ccccc1", style="plain", canvas_size=800, padding=20.0)
# Pass output_path="out.svg" to also write a file.
```

## Drawing styles

| Style | Description |
|-------|-------------|
| `plain` | Colored disks; carbon atoms unlabeled |
| `names_hetero` | Element symbols on heteroatoms |
| `names_all` | Symbols on every atom |
| `stroke` | Filled disks with bond-colored outline |

Use `--style` on the CLI or the `style=` argument in `render_smiles_to_svg`.

### Examples (default style)

```bash
smiles2svg -s 'Cn1cnc2n(C)c(=O)n(C)c(=O)c12'
```

![Caffeine — default](figures/caffeine_default.svg)

### Stroke, labels, color, implicit H on donors

<p align="middle">
  <img src="figures/caffeine_stroke.svg" width="24%" alt="stroke" />
  <img src="figures/caffeine_names_hetero.svg" width="24%" alt="names hetero" />
  <img src="figures/caffeine_names_all.svg" width="24%" alt="names all" />
  <img src="figures/caffeine_color.svg" width="24%" alt="custom color" />
</p>

```bash
smiles2svg -s 'Cn1cnc2n(C)c(=O)n(C)c(=O)c12' --style stroke
smiles2svg -s 'Cn1cnc2n(C)c(=O)n(C)c(=O)c12' --style names_hetero
smiles2svg -s 'Cn1cnc2n(C)c(=O)n(C)c(=O)c12' --style names_all
smiles2svg -s 'Cn1cnc2n(C)c(=O)n(C)c(=O)c12' --style names_hetero --color '#468BFF'
```

### Implicit hydrogens on H-bond donors

```bash
smiles2svg -s 'Cn1cnc2n(C)c(=O)n(C)c(=O)c12' --add-hydrogens
```

![Add hydrogens](figures/add_hydrogens.svg)

## Documentation

The manual (**Furo** theme) covers getting started (install, CLI), user guide (styles, canvas, Python API), examples gallery, citation, and API reference — hosted on [Read the Docs](https://smiles2svg.readthedocs.io).

Build locally:

```bash
pip install -e ".[docs]"
cd docs && make html
```

Regenerate embedded SVG figures (optional):

```bash
python docs/scripts/generate_doc_figures.py
```

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md).

## License

MIT — see [LICENSE](LICENSE).

## Roadmap (from original project)

- [x] Custom font
- [x] PNG export
- [x] Custom basename for output files
- [x] Stroke style for all atoms
- [x] Optional implicit H on typical donors
- [ ] Improve text centering for labels
- [ ] Bond color gradients between atoms
- [ ] Lone pairs
