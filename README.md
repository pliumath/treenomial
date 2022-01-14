
# treenomial 1.2.5

## Overview

The package **treenomial** is an application of polynomials that
uniquely describe trees. It provides tools for tree analysis and
comparison based on polynomials. The core functions are:

  - **`treeToPoly()`**: convert rooted unlabeled binary trees to tree
    distinguishing polynomials described with coefficient matrices

  - **`polyToDistMat()`**: construct a distance matrix from multiple
    coefficient matrices using a distance measure

For the mathematical description of the tree defining polynomial see:

[Liu, Pengyu. “A tree distinguishing polynomial.” arXiv preprint
arXiv:1904.03332 (2019).](https://arxiv.org/abs/1904.03332)

## Installation

    library(devtools)
    install_github("pliumath/treenomial")
