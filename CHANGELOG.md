# Changelog

## 0.2.0 — Unreleased

### Changed

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
