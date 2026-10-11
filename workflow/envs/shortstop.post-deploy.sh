#!/usr/bin/env bash
# Run by Snakemake inside the new env, after shortstop.yaml is installed.
# Bump the commit to update ShortStop (146d40c = release v1.1.0).
#
# --use-pep517 builds in an isolated env with a fresh setuptools. A legacy
# `setup.py egg_info` inside this env pairs setuptools 80 with an older
# `packaging` and dies in canonicalize_version(strip_trailing_zero=...).
set -euo pipefail
SHORTSTOP_COMMIT=146d40c11d2db1458001b912d11d949684a1f206
python -c "import setuptools, packaging; print('env: setuptools', setuptools.__version__, 'packaging', packaging.__version__)"
pip install --no-deps --use-pep517 "git+https://github.com/brendan-miller-salk/ShortStop.git@${SHORTSTOP_COMMIT}"
python -c "import shortstop, umap, xgboost; print('shortstop import OK')"
