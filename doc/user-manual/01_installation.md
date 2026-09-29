# Chapter 1 -- Installation

## Table of contents

- [Overview](#overview)

## Overview

Spinniped requires Python 3.13 or newer, NumPy, SciPy, and Matplotlib. These
runtime dependencies are installed automatically with the package.

Until a release is available on PyPI, install a local clone in editable mode:

```bash
git clone https://github.com/davidepisu9723/spinniped.git
cd spinniped
python -m pip install -e .
```

Install the test dependency with:

```bash
python -m pip install -e ".[test]"
```
