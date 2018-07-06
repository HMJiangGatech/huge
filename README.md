# huge
High-Dimensional Undirected Graph Estimation

## Installation

## Reference
[1] [High-dimensional Undirected
Graph Estimation](https://cran.r-project.org/web/packages/huge/vignettes/vignette.pdf)  
[2] [Quanquan Gu, Yuan Cao, et al. Local and Global Inference for High Dimensional Nonparanormal Graphical Models](https://arxiv.org/abs/1502.02347)  
[3] [Picasso: A Sparse Learning Library for High Dimensional Data Analysis in R and Python](https://cran.r-project.org/web/packages/picasso/vignettes/vignette.pdf)

## GSOC TODO List

### Code Refactoring and Benchmarking
- [x] Use `roxygen2` to manage Documents.
- [ ] Benchmark the code.
- [ ] Use `RcppEigen` to accelerate the code.
- [ ] Use `OpenMP` to enable multi-processing.

### Solver Upgrade
- [ ] Implement Active Set Newton solver.

### Inference Module
- [ ] Add inference module

### Polish and submit
- [ ] Polish `readme.md`
- [ ] Polish R Documents
- [ ] Polish `vignette`
- [ ] submit to CRAN
