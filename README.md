# Computing invariants of persistence modules

We provide a code implementation for computing 
1. Interval rank invariants of persistence modules under total or source-sink compression systems. Each interval rank (also called compression multiplicity **[2]**) is computed with the formula provided in **[3]** via rank computation of structure linear maps.
2. Interval replacement of persistence modules under total or source-sink compression systems. Each signed interval multiplicity is computed using Möbius inversion of interval ranks. 
3. Interval multiplicities of persistence modules, that is, the multiplicities of interval modules as direct summands in persistence modules. Each interval multiplicity is computed with the formula provided in **[4]** via rank computation of structure linear maps.

This code can also find the maximal interval-decomposable summand of a given persistence module. Thus, one can determine the interval decomposability by this code.

We also provide ways of visualization by using standard persistence diagrams for 1D persistence modules (representations) and connected persistence diagrams (proposed in **[5]**) for ladder representations.

## Background

In topological data analysis, the usual pipeline is the following:
1. Geometric data
2. Filtered simplicial complexes (obtained with a choice of filtration)
3. Persistence module M (obtained with a choice of homology functor)
4. **Algebraic signature of M (obtained with a choice of algebraic descriptor)**
5. Vector that can be used in machine learning tasks (obtained with a choice of vectorization)

Tasks 1 to 3 should be straightforward to implement thanks to existing code on other repositories. This repository focuses on **task 4**. Given a persistence module, the code computes its interval replacement, used as an algebraic descriptor. This choice of descriptor is motivated by **[1]**: the interval replacement preserves the so-called interval ranks, which generalize the rank invariant property of other algebraic descriptors. Task 5 is still an open problem for the interval replacement.

## Overview

In this repository, you will find:
- three **tutorial** notebooks, in which we explain how to compute interval replacements, interval ranks, interval multiplicities, and other features;
- a **utils.py** file containing all the necessary code implementations used in the notebook ;
-  and also a **display.py** file containing some useful functions to visualize representations in 1D or 2D settings.

## Key features

In the notebook, you will find how to:
- **instantiate a d-dimensional grid** which is the quiver considered here. This is done within the class `Representation` ; 
- **define a representation (persistence module)** by adding vector spaces and linear maps to the quiver ;
- **define intervals** of the quiver with the class `Interval`. By default, intervals are defined by a list of sources and a list of sinks. One can access all points within the interval using `int_hull`. Conversely, given a list of points forming a connected and convex set, one can instantiate an `Interval` object by using `get_src_snk` ;
- **obtain the list of all intervals** thanks to `list_int` ;
- **compute the interval rank** of a given interval. It is computed with the formula from **[3]** via rank computation of linear maps ;
> [!NOTE]
> Rank computations are done in $\mathbb{R}$. You can adjust the rank computation algorithm to handle rank computation in finite fields for example. 
- **compute the interval signed multiplicity** of a given interval via Möbius inversion, by computing the cover of the interval. Signed multiplicities yield the interval replacement of the persistence module.

Additionally, we provide some **visualization** features for the quiver and its intervals.

By default "total (tot)" compression system is used. To use "source-sink (ss)" instead, add a `compression` argument: `R.int_replacement(interval, compression='ss')`.

## Installation

This implementation is built from scratch and does not depend on any external Python libraries except SymPy. You can install it with `pip` by running the following command in a terminal:
```
pip install sympy
```

To use the code, download **all the files** from this repository and place them in the **same folder**.
You can then open and run the notebook **tutorial.ipynb**, which demonstrates how the code works on a specific example of a quiver representation.

From there, you can easily modify the code in the notebook to compute the interval replacements for any representation you are interested in.

## Usage

This code is distributed in the hope that it will be useful. It might be integrated into a topological data analysis pipeline to provide algebraic descriptors from data directly. We are developing the code and will add more new features in the future.

## References

**[1]**: Asashiba, H., Buchet, M., Escolar, E. G., Nakashima, K., & Yoshiwaki, M. *On interval decomposability of 2D persistence modules*, Computational Geometry, Volumes 105–106, 2022, 101879, ISSN 0925-7721, [https://doi.org/10.1016/j.comgeo.2022.101879](https://doi.org/10.1016/j.comgeo.2022.101879).

**[2]**: Asashiba, H., Escolar, E. G., Nakashima, K., & Yoshiwaki, M. *On Approximation of 2D Persistence Modules by Interval-decomposables*. Journal of Computational Algebra, Volumes 6–7, 2023, 100007, ISSN 2772-8277, [https://doi.org/10.1016/j.jaca.2023.100007](https://doi.org/10.1016/j.jaca.2023.100007).

**[3]**: Asashiba, H., Gauthier, E., & Liu, E. *Interval Replacements of Persistence Modules*. arXiv preprint [arXiv:2403.08308](https://arxiv.org/abs/2403.08308) (2024).

**[4]**: Asashiba, H., & Liu, E. *Interval Multiplicities of Persistence Modules*. arXiv preprint [arXiv:2411.11594](https://doi.org/10.48550/arXiv.2411.11594) (2024).

**[5]**: Hiraoka, Y., Nakashima, K., Obayashi, I., & Xu, C. *Refinement of interval approximations for fully commutative quivers*. Japan Journal of Industrial and Applied Mathematics, Volumes 42, Issue 4, 1309-1361 (2025), [https://doi.org/10.1007/s13160-025-00739-w](https://doi.org/10.1007/s13160-025-00739-w).

**[6]**: Kim, W., & Mémoli, F. *Generalized persistence diagrams for persistence modules over posets*. J Appl. and Comput. Topology 5, 533–581 (2021). [https://doi.org/10.1007/s41468-021-00075-1](https://doi.org/10.1007/s41468-021-00075-1).