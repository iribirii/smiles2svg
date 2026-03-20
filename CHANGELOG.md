# Changelog

## 0.2.0 — Unreleased

### Changed

- **Python:** Requires **3.10+** (was 3.9+).
- **Docs:** Sphinx manual uses the **Furo** theme with sections for getting started, user guide, examples gallery (generated SVGs), citation, and API. Script `docs/scripts/generate_doc_figures.py` rebuilds gallery images.
- **Rendering:** SVG is built in **pixel coordinates** (like `xyzrender`): `viewBox="0 0 W H"` with explicit `width`/`height`, `font-size` in **px**. Bond stroke and label size are **Å × scale** (fixed proportion to RDKit coordinates), not tied to a reference span — canvas- and span-independent proportions. New `layout` module; CLI/API options `--canvas-size` / `--padding` and `canvas_size` / `padding` kwargs.
- **Packaging:** Switched to `pyproject.toml`, `src/` layout, and entry point `smiles2svg`.
- **Structure:** Split modules into `render`, `drawing`, `molecule`, `elements`, and `cli`.
- **API:** Added `render_smiles_to_svg()` for library use; CLI is a thin wrapper.
- **CLI:** `--smiles` / `--smiles-file` are explicit; `--smiles-file` is preferred over magic `'none'` strings.
- **CLI:** Custom `--name` applies only when a single SMILES is passed (avoids overwriting output when batching).

### Fixed

- Invalid SMILES now raise `ValueError` instead of failing later.
- Element color table matches the original implementation byte-for-byte.

## 0.1.1 and earlier

- Original script-style layout (`setup.py`, flat package without `__init__.py`).
