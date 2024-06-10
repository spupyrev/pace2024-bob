# Bob
This repository contains the source code of a solver for the [PACE 2024](https://pacechallenge.org/2024/) challenge.

## Description

Coming

## Requirements
Your machine needs to have a C++ compiler with C++17 support. We use GNU Make to build the executables.

The solver has been developed and thoroughly tested with LLVM ([release_17.0.4](https://releases.llvm.org/)). 
Use the corresponding versions of clang/lld for reproducible results.

## Building & Running

To build the executable, run one of the following commands at the root of the repository:
```
make heuristic -j
make exact -j
make cutwidth -j
```
This creates an executable named `pace_*` and can be run as follows:

| Track      | Command |
| -------------- | ------- |
| heuristic  | `./pace_heuristic < path/to/input_file.gr` |
| exact      | `./pace_exact < path/to/input_file.gr` |
| cutwidth   | `./pace_cutwidth < path/to/input_file.gr` |
| dev (`make -j`)       | `./pace < path/to/input_file.gr` or `./pace -verbose=1 -help` for supported options |

> [!IMPORTANT]
> Release **pace-2024** contains statically built binaries (for Linux) for the three tracks: `pace_heuristic`, `pace_exact`, and `pace_cutwidth`
