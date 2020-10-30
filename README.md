# MMSP Nucleation Benchmark

This repository contains an implementation of the [PFHub Nucleation
Benchmark](https://pages.nist.gov/pfhub/benchmarks/benchmark8.ipynb/)
using [MMSP](https://github.com/mesoscale/mmsp).

## 1-explicit

This folder implements Problem 8.1: Explicit Nucleation.

## 2-initial-nucleation

This folder implements Problem 8.2: Random Nucleation at *t*=0. It uses
the [CImg](http://cimg.eu) library to determine the number of disjoint
particles, after [this answer](https://stackoverflow.com/a/40898020) on
StackOverflow.

## 3-initial-nucleation

This folder implements Problem 8.3: Random Nucleation at Random Times.
It uses the [CImg](http://cimg.eu) library to determine the number of
disjoint particles, after [this answer](
https://stackoverflow.com/a/40898020) on StackOverflow.

## 4-athermal

This folder implements Problem 8.4: Athermal Nucleation.

## 7-initial-nucleation

This folder modifies Problem 8.3 to add a single nucleus at *t*=0.

---

CImg is distributed under the [CeCILL](
https://opensource.org/licenses/CECILL-2.1) open-source software
license. The remainder of this work was conducted by an employee of the
Federal Government during the course of their normal duties, and is not
subject to copyright protection within the United States of America.
