
# treenomial 1.2.5

## Overview

The package **treenomial** is an application of polynomials that
uniquely describe trees and a tree lattice that serves as a coordinate system for rooted binary trees. It provides tools for tree analysis and
comparison based on polynomials and the tree lattice. The core functions are:

  - **`treeToPoly()`**: convert rooted unlabeled binary trees to tree
    distinguishing polynomials described with coefficient matrices

  - **`treeToLattice()`**: construct rooted binary trees with branch lengths 
    to the lattice representation

For the mathematical description of the tree defining polynomial see:

[Liu, Pengyu, A tree distinguishing polynomial, *Discrete Applied Mathematics* 288 (2021).](https://doi.org/10.1016/j.dam.2020.08.019)

## Installation

    library(devtools)
    install_github("pliumath/treenomial")
