# Bob
This repository contains the source code of a solver for the [PACE 2024](https://pacechallenge.org/2024/) challenge.

## Description

Coming

## Requirements
Your machine needs to have a C++ compiler with C++17 support (Clang is recommended but GCC works too). 
We use GNU Make to build the executables.

## Building & Running

To build the executable, run the following commands at the root of the repository:
```
make r -j
```
This creates an executable named `pace_release` and can be run as follows:

| Track      | Command |
| -------------- | ------- |
| heuristic (*default*)  | `./pace_release < path/to/input_file.gr` |
| exact      | `./pace_release -confidence=30 -time-limit=600 < path/to/input_file.gr` |
| cutwidth   | `./pace_release -time-limit=30 < path/to/input_file.gr` |
