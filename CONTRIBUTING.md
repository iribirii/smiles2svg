# Contributing

## Development setup

```bash
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -e ".[dev,docs]"
```

## Tests

```bash
pytest
```

Tests assume RDKit and other runtime dependencies are installed (same as a normal install).

## Documentation

The manual uses **Sphinx** with the **Furo** theme (`pip install -e ".[docs]"`).

```bash
cd docs && make html
```

Open `docs/_build/html/index.html`. After changing how diagrams render, regenerate the SVG figures used in the docs:

```bash
python docs/scripts/generate_doc_figures.py
```

## Style

- Match existing formatting and naming in the files you touch.
- Prefer small, focused changes tied to a single issue or feature.

## Pull requests

1. Describe the change and link to an issue if one exists.
2. Ensure tests pass and documentation builds.
3. Keep commits scoped and messages clear.
