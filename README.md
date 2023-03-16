# fluxome

<p align="center">
    <em>integrated probabilistic modeling of stochastic gene regulatory networks and cell state trajectories</em>
</p>

[![build](https://github.com/biotheorylab/fluxome/workflows/Build/badge.svg)](https://github.com/biotheorylab/fluxome/actions)
[![codecov](https://codecov.io/gh/biotheorylab/fluxome/branch/main/graph/badge.svg)](https://codecov.io/gh/biotheorylab/fluxome)
[![PyPI version](https://badge.fury.io/py/fluxome.svg)](https://badge.fury.io/py/fluxome)

---

**Documentation**: <a href="https://biotheorylab.github.io/fluxome/" target="_blank">https://biotheorylab.github.io/fluxome/</a>

**Source Code**: <a href="https://github.com/biotheorylab/fluxome" target="_blank">https://github.com/biotheorylab/fluxome</a>

---

**Table of Contents**

- [fluxome](#fluxome)
  - [Installation](#installation)
  - [Development](#development)
    - [Specification](#specification)
    - [Setup environment](#setup-environment)
    - [Run unit tests](#run-unit-tests)
    - [Format the code](#format-the-code)
    - [Publish a new version](#publish-a-new-version)
  - [Serve the documentation](#serve-the-documentation)
  - [License](#license)

## Installation

```console
python -m pip install "fluxome @ git+https://github.com/biotheorylab/fluxome.git@main"
```

## Development

### Specification

See [SPECIFICATION.md](./SPECIFICATION.md).

### Setup environment

We use [Hatch](https://hatch.pypa.io/latest/install/) to manage the development environment and production build. Ensure it's installed on your system. It is often convenient to do this with [pipx](https://pypa.github.io/pipx/installation/).

You can print the environments and script commands supported by the current content of [pyproject.toml](./pyproject.toml) with:

```bash
hatch env show
```

### Run unit tests

You can run all the tests with:

```bash
hatch run test
```

### Format the code

Execute the following command to apply linting and check typing:

```bash
hatch run lint
```

### Publish a new version

You can bump the version, create a commit and associated tag with one command:

```bash
hatch version patch
```

```bash
hatch version minor
```

```bash
hatch version major
```

Your default git text editor will open so you can add information about the release.

When you push the tag on GitHub, the workflow will automatically publish it on PyPi and a GitHub release will be created as draft.

## Serve the documentation

You can serve the Mkdocs documentation with:

```bash
hatch run docs-serve
```

This will automatically watch for changes in your code.

## License

This project is licensed under the terms of the MIT license.
