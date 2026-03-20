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

Build the Sphinx site locally:

```bash
cd docs
make html
```

Open `docs/_build/html/index.html` in a browser.

## Style

- Match existing formatting and naming in the files you touch.
- Prefer small, focused changes tied to a single issue or feature.

## Pull requests

1. Describe the change and link to an issue if one exists.
2. Ensure tests pass and documentation builds.
3. Keep commits scoped and messages clear.
