---
jupyter:
  celltoolbar: Slideshow
  jupytext:
    cell_metadata_json: true
    formats: ipynb,md,py:percent
    notebook_metadata_filter: all
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.1
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
  language_info:
    codemirror_mode:
      name: ipython
      version: 3
    file_extension: .py
    mimetype: text/x-python
    name: python
    nbconvert_exporter: python
    pygments_lexer: ipython3
    version: 3.10.7
  rise:
    scroll: true
    theme: black
  toc-autonumbering: true
  toc-showcode: false
  toc-showmarkdowntxt: false
  widgets:
    application/vnd.jupyter.widget-state+json:
      state: {}
      version_major: 2
      version_minor: 0
---

# Wright-Fisher model

```python tags=[]
from typing import List

import numpy as np
```

```python tags=[]
def wright_fisher_fwd(
    N: int = 100,
    T: int = 50,
    nSim: int = 20,
    theta_f: List[float] = [0.5],
    theta_h: List[float] = [0.05],
    theta_z0: List[float] = [0.1],
    seed: float = 100,
    verbose: bool = False,
) -> (np.ndarray, np.ndarray, np.ndarray):
    """Forward simulation of the Wright-Fisher model.

    Args:
        N (int, optional): size of population. Defaults to 100.
        T (int, optional): time-points. Defaults to 50.
        nSim (int, optional): number of forward simulations. Defaults to 20.
        theta_f (List[float], optional): log relative fitness(es) of given variant(s). Defaults to a single variant with 0.5.
        theta_h (List[float], optional): mutation rate of variant(s). Defaults to single variant with mutation rate 0.05.
        theta_z0 (List[float], optional): initial probability of variant(s). Defaults to single variant with probability 0.1.
        seed (float, optional): random seed. Defaults to 100.
        verbose (bool, optional): increase printing verbosity. Defaults to False.

    Returns:
        Zs, Pis, log_Ps ((numpy.ndarray, numpy.ndarray, numpy.ndarray)): _description_

    Examples:
        To run with default paramaters

            >>> Zs, Pis, log_Ps = fisher_wright_fwd()

    Todo:
        * Upgrade nested List outputs to numpy.ndarray
    """

    Zs = [None] * nSim
    Pis = [None] * nSim  # np.empty((nSim))
    log_Ps = np.zeros((nSim))
    Dz = len(theta_f)

    np.random.seed(seed)

    for cSim in range(nSim):

        if verbose:
            print(f"\nsimulation: {cSim}")

        Z = np.zeros((T, N, Dz))
        Pi = np.zeros((T, N))
        for i in range(Dz):
            Z[0, :, i] = np.random.rand(1, N) < theta_z0[i]
            log_Ps[cSim] = np.sum(
                Z[0, :, i] * np.log(theta_z0[i])
                + (1 - Z[0, :, i]) * np.log(1 - theta_z0[i])
            )

        for t in range(1, T):
            fs = np.zeros((N))
            for i in range(Dz):
                fs = fs + Z[t - 1, :, i] * theta_f[i]

            fs = np.exp(fs)
            fs = fs / np.sum(fs)
            fs_ = np.cumsum(fs)

            for n in range(N):
                idx = np.argwhere(np.random.rand() <= fs_).flatten()
                idx = idx[0]
                Pi[t, n] = idx
                log_Ps[cSim] = log_Ps[cSim] + np.log(fs[idx])
                for i in range(Dz):
                    if np.random.rand() >= theta_h[i]:
                        Z[t, n, i] = Z[t - 1, idx, i]
                        log_Ps[cSim] = log_Ps[cSim] + np.log(1 - theta_h[i])
                    else:
                        Z[t, n, i] = 1 - Z[t - 1, idx, i]
                        log_Ps[cSim] = log_Ps[cSim] + np.log(theta_h[i])

        Zs[cSim] = Z
        Pis[cSim] = Pi

        if verbose:
            # [print(x.shape) for x in (Z, Pi, log_Ps)]
            # [print(x.shape) for x in (fs,)]
            # [print(np.array(x, dtype=object).shape) for x in (Zs, Pis)]
            for key, value in dict(
                [
                    ("Z", Z),
                    ("Pi", Pi),
                    ("log_Ps", log_Ps),
                    ("fs", fs),
                    ("Zs", Zs),
                    ("Pis", Pis),
                ]
            ).items():
                print(f"size {key}: {np.array(value, dtype=object).shape}")

    return Zs, Pis, log_Ps
```

```python tags=[]
# fwd simulations
nSim1 = 3
# bkw simulations
nSim2 = 5

# time-points
T = 4

# size of population
N = 5

# log relative fitness of variant
theta_f = [0.5]
# mutation rate
theta_h = [0.05]
# initial probability of variant
theta_z0 = [0.1]

# mixing param for bkw simulation
alpha = 0.85 * np.ones((1, T))
# verbosity
verbose = True
# for smoothing proposal dist
ep = 0

seed = 100
```

```python tags=[]
Zs, Pis, log_Ps = wright_fisher_fwd(
    N, T, nSim1, theta_f, theta_h, theta_z0, seed, verbose
);
```

```python
type(Zs)
```

```python tags=[]
assert np.array(Zs).shape == (nSim1, T, N, len(theta_f))
```

```python tags=[]
assert np.array(Pis).shape == (nSim1, T, N)
```

```python tags=[]
assert np.array(log_Ps).shape == (nSim1,)
```

```python tags=[]
for key, value in dict([("Zs", Zs), ("Pis", Pis), ("log_Ps", log_Ps)]).items():
    print(f"size {key}: {np.array(value).shape}")
```

```python

```