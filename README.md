<h1 align="center">Huge</h1>

[![](https://cranlogs.r-pkg.org/badges/huge)](https://cran.r-project.org/package=huge)[![](https://cranlogs.r-pkg.org/badges/grand-total/huge)](https://cran.r-project.org/package=huge)

<h4 align="center">R Package for High-Dimensional Undirected Graph Estimation and Inference</h4>

___Huge___ (Huge-Dimensional Undirected Graph Estimation) implements the algorithm of estimating the parameters of a Gaussian distribution in such a way that the resulting undirected graphical model is sparse. The core algorithm is implemented in C++ with RcppEigen support for portable high performance linear algebra. It also implements a unified framework to quantify local and global inferential uncertainty for high dimensional graphical models. In particular, we consider the problems of testing the presence of a single edge. Runtime profiling is documented in the [__Performance__](#performance) section.


## Installation

### Prerequisites

Huge uses OpenMP to enables faster matrix multiplication. So, to use huge, you must correctly enables OpenMP for the compiler.

For Windows and Linux users, newest version of GCC has fully support of OpenMP.

But for MAC OS users, things are a little tricky since the default llvm on MAC OS does not support OpenMP. But the solution is easy. You can simply install llvm with full OpenMP support and direct R using this version of llvm.

First, install llvm with OpenMP support by typing

```
brew install llvm
```

Then append the following lines into `~/.R/Makevars` to enable llvm with OpenMP support to be the compiler for R packages.

```
CC = /usr/local/bin/clang-omp
CXX = /usr/local/bin/clang-omp++
CXX98 = /usr/local/bin/clang-omp++
CXX11 = /usr/local/bin/clang-omp++
CXX14 = /usr/local/bin/clang-omp++
CXX17 = /usr/local/bin/clang-omp++
OBJC = /usr/local/bin/clang-omp
OBJCXX = /usr/local/bin/clang-omp++
```

### Installing from GitHub

First, you need to install the devtools package. You can do this from CRAN. Invoke R and then type

```
install.packages(devtools)
```

Then load the devtools package and install huge

```
library(devtools)
install_github("HMJiangGatech/huge")
library(huge)
```

*Windows User:*  If you encounter a Rtools version issue: 1. make sure you install the latest [Rtools](https://cran.r-project.org/bin/windows/Rtools/); 2. try the following code
```R
assignInNamespace("version_info", c(devtools:::version_info, list("3.5" = list(version_min = "3.3.0", version_max = "99.99.99", path = "bin"))), "devtools")
```

### Install from CRAN

Ideally you can just install and enable huge using with the help of CRAN on an R console.

```
install.packages("huge")
library(huge)
```

## Examples

```R
#generate data  
L = huge.generator(n = 50, d = 12, graph = "hub", g = 4)

#graph path estimation using glasso  
est = huge(L$data, method = "glasso")
plot(est)

#inference of Gaussian graphical model at 0.05 significance level  
T = est$icov[[10]]  
inf = huge.inference(L$data, T, L$theta)
print(inf$error) # print out type-I error
```

## Experiments
For detailed implementation of the experiments, please refer to `benchmark/benchmark.R`

### Graph Estimation

We compared our package on hub graph with (n=200,d=200) with other packages, namely, QUIC and clime.
Huge significantly outperforms clime, QUIC and original huge in timing performance. We also calculated the likelihood for estimation.

<center>
<table>
  <thead>
    <tr>
      <th></th>
      <th>CPU Times(s)</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><center>Huge glasso</center></td>
      <td><center>1.12</center></td>
    </tr>
    <tr>
      <td><center>Huge tiger</center></td>
      <td><center>1.88</center></td>
    <tr>
      <td><center>Huge v1.2.7</center></td>
      <td><center>1.80</center></td>
    </tr>
	<tr>
	  <td><center>QUIC</center></td>
      <td><center>7.50</center></td>
	</tr>
	<tr>
	  <td><center>Clime</center></td>
      <td><center>416.77</center></td>
	</tr>
  </tbody>
</table>
</center>

<center>
<table>
  <thead>
    <tr>
      <th></th>
      <th>Object value</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><center>Huge glasso</center></td>
      <td><center>-125.96</center></td>
    </tr>
    <tr>
      <td><center>Huge tiger</center></td>
      <td><center>-125.47</center></td>
	<tr>
	  <td><center>QUIC</center></td>
      <td><center>-90.58</center></td>
	</tr>
	<tr>
	  <td><center>Clime</center></td>
      <td><center>-136.96</center></td>
	</tr>
  </tbody>
</table>
</center>

### Graph Inference
When using the Gaussian graphical model, huge controls the type I error well.

<center>
<table>
    <tr>
	  <td></td>
	  <td colspan = "2"><center>band</center></td>
	  <td colspan = "2"><center>hub</center></td>
	  <td colspan = "2"><center>scale-free</center></td>
	<tr>
      <td><center>significance level</center></td>
      <td><center>0.05</center></td>
      <td><center>0.10</center></td>
	  <td><center>0.05</center></td>
      <td><center>0.10</center></td>
      <td><center>0.05</center></td>
      <td><center>0.10</center></td>
    </tr>
    <tr>
      <td><center>type I error</center></td>
	  <td><center>0.0175</center></td>
	  <td><center>0.0391</center></td>
      <td><center>0.0347</center></td>
	  <td><center>0.0669</center></td>
	  <td><center>0.0485</center></td>
      <td><center>0.0854</center></td>
	</tr>
  </tbody>
</table>
</center>

## References
[1] [T. Zhao and H. Liu, The huge Package for High-dimensional Undirected Graph Estimation in R, 2012](https://cran.r-project.org/web/packages/huge/vignettes/vignette.pdf)  
[2] [Xingguo Li, Jason Ge, Haoming Jiang, Mingyi Hong, Mengdi Wang, and Tuo Zhao, Boosting Pathwise Coordinate Optimization: Sequential Screening and Proximal Subsampled Newton Subroutine, 2016](https://www2.isye.gatech.edu/~tzhao80/)  
[3] [Quanquan Gu, Yuan Cao, et al. Local and Global Inference for High Dimensional Nonparanormal Graphical Models](https://arxiv.org/abs/1502.02347)  
[4] [Conﬁdence intervals for high-dimensional inverse covariance estimation](https://projecteuclid.org/download/pdfview_1/euclid.ejs/1433195859)  
[5] D. Witten and J. Friedman, New insights and faster computations for the graphical lasso,2011  
[6] N. Meinshausen and P. Buhlmann, High-dimensional Graphs and Variable Selection with the Lasso, 2006
