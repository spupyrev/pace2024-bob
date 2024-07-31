# Bob
This repository contains the source code of a solver for the [PACE 2024](https://pacechallenge.org/2024/) challenge.
See a [short description](docs/pace24bob.pdf) of the solver.

Bob came third&#129353;on the exact track, second&#129352;on the heuristic track, and got the maximum score on the parameterized track.
    
## Requirements
Your machine needs to have a C++ compiler with C++17 support. We use GNU Make to build the executables.

The solver has been developed and thoroughly tested with LLVM ([release_17.0.4](https://releases.llvm.org/)). 
Use the corresponding versions of clang/lld for reproducible results.

## Building & Running

To build the executable, run one of the following commands at the root of the repository:
```
make heuristic
make exact
make cutwidth
```
This creates an executable named `pace_*` and can be run as follows:

| Track      | Command |
| -------------- | ------- |
| heuristic  | `./pace_heuristic < path/to/input_file.gr` |
| exact      | `./pace_exact < path/to/input_file.gr` |
| cutwidth   | `./pace_cutwidth < path/to/input_file.gr` |
| dev (`make`)       | `./pace < path/to/input_file.gr` or `./pace -verbose=1 -help` for supported options |
| lite (`make lite`)       | `./pace_lite < path/to/input_file.gr` |

> [!IMPORTANT]
> Release **pace-2024** contains statically built binaries (for Linux) for the three tracks: `pace_heuristic`, `pace_exact`, and `pace_cutwidth`

The lite version of the solver, `pace_lite` (not part of the competition), is able to process very large instances within seconds 
at the cost of slightly worsened results. For example, *every* instance from the public dataset (containing up to 100K vertices) is processed 
within 10 seconds producing less than extra 0.5% of crossings than the known optimum.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11551133.svg)](https://doi.org/10.5281/zenodo.11551133)
