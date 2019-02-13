## Test environments
* local OS X install, R 3.1.2
* ubuntu 16.04, R 3.5.2
* windows 10, R 3.5.2

## R CMD check results

I got one NOTE when checking under ubuntu environment:
```
* checking installed package size ... NOTE
  installed size is 10.9Mb
  sub-directories of 1Mb or more:
    libs   9.2Mb
```
It seems that on LINUX architectures, the CHECK returns one NOTE because the libs subdirectory is then above the 1MB threshold. However, it seems that this NOTE only appears under LINUX, but not under Windows or OSX.
My understanding is that this inflation of the libs subdirectory is due to the use of Rcpp. Indeed, some functions of the package have been written in C++ using Rcpp. They are needed to perform [what the package do]. Without the speed up gained from those C++ functions, this package would become impractical.


## Downstream dependencies
