## Test environments
* local OS X install, R 3.1.2
* ubuntu 16.04, R 3.5.2
* windows 10, R 3.5.2

## R CMD check results

One WARNING under ubuntu environment:
```
* checking compilation flags used ... WARNING
  Compilation used the following non-portable flag(s):
    ‘-Wdate-time’ ‘-Werror=format-security’ ‘-Wformat’
```
This warning only appears under ubuntu environment. This is because the default compiler flag setting is specical on the system. Problem solved by https://stackoverflow.com/questions/50658198/cran-submission-r-cmd-check-warning-compilation-flags-used.

One NOTE under ubuntu environment:
```

* checking installed package size ... NOTE
  installed size is 10.9Mb
  sub-directories of 1Mb or more:
    libs   9.2Mb
```
It seems that on LINUX architectures, the CHECK returns one NOTE because the libs subdirectory is then above the 1MB threshold. However, it seems that this NOTE only appears under LINUX, but not under Windows or OSX. My understanding is that this inflation of the libs subdirectory is due to the use of Rcpp.


## Downstream dependencies
I have also run R CMD check on downstream dependencies of huge.
All packages that I could install passed, except for `netgwas`:
  One check is failed: testthat/test_buildMap.R
  I manually verified the difference between the output of the old huge and the new huge. The only difference is the ordering of the result.
  There is only one call of our `huge` function. I also check the output of that step. I only discovered subtle numerical difference. The relative error is about $10^{-10}$.
  
