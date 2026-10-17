# Contributing to PPDyn

Issues and pull requests are welcome on GitHub.

## Development setup

```bash
git clone https://github.com/sayanadhikari/PPDyn.git
cd PPDyn
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
pytest -q
```

## Pull requests

1. Branch from `main`. New physics belongs on `main` behind an input switch, not on a long-lived feature branch.
2. Add tests for new behaviour; physics features should be checked against an analytic or published result.
3. Keep `pytest -q` passing.

## Adding a feature

* **New input keys** go in `src/ppdyn/schema.py` with a default, type, physical dimension and one-line description. Validation, SI conversion and `docs/source/input_reference.md` (regenerate with `ppdyn docs -o docs/source/input_reference.md`) all follow from that entry.
* **Kernels** live in `src/ppdyn/kernels.py`. Parallel loops must be race-free: each particle writes only its own entries, as in `_pair_sorted` and `collide`.
* **External forces** are packed into the parameter vector in `Simulation._setup_forces` and evaluated in `kernels.external_forces`.
* **Examples** go in `src/ppdyn/examples/<name>/input.toml`, with the first comment line used as the description in `ppdyn init`.

## Coding style

PEP 8, 4-space indentation, lines up to about 110 characters.

## Issues

Please include the input file, the PPDyn version (`ppdyn --version`) and the output of `ppdyn check input.toml`.

## License

By contributing to PPDyn you agree that your contributions are licensed under the MIT License.
