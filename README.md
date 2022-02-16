# `XtalsPyTools.jl`

| **Build Status** | **Test Coverage** |
|:---:|:---:|
| [![Build](https://github.com/SimonEnsemble/XtalsPyTools.jl/actions/workflows/ci_testing.yml/badge.svg)](https://github.com/SimonEnsemble/XtalsPyTools.jl/actions/workflows/ci_testing.yml) | [![codecov](https://codecov.io/gh/SimonEnsemble/XtalsPyTools.jl/branch/master/graph/badge.svg?token=JD0C2Y99H1)](https://codecov.io/gh/SimonEnsemble/XtalsPyTools.jl) |

A repository of Python-dependent tools for working with [`Xtals.jl`](https://github.com/SimonEnsemble/Xtals.jl)

`Xtals` is a pure-Julia package for reading/writing/editing crystal structures.
Some higher-level features were previously dependent on Python packages; those have been moved here.

To obtain the package, simply add it via the REPL package manager:

`pkg> add XtalsPyTools`

The methods in this package require `scipy` and/or `pymatgen` to be installed.
To configure these dependencies automatically in a Julia-managed environment, use:

`julia quick_setup.jl`

See the [`Xtals.jl` docs](https://SimonEnsemble.github.io/Xtals.jl/stable) for further information.
