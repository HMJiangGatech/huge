## Test environments
* Mac OS X, R 3.5.2
* ubuntu 16.04, R 3.5.2
* windows 10, R 3.5.2

## R CMD check results

One NOTE under ubuntu environment:
```

* checking installed package size ... NOTE
  installed size is 10.9Mb
  sub-directories of 1Mb or more:
    libs   9.2Mb
```
It seems that on LINUX architectures, the CHECK returns one NOTE because the libs subdirectory is then above the 1MB threshold. However, it seems that this NOTE only appears under LINUX, but not under Windows or OSX. My understanding is that this inflation of the libs subdirectory is due to the use of Rcpp.


## Downstream dependencies
I have also run R CMD check on downstream dependencies of huge by https://github.com/r-lib/revdepcheck/.
All 21 packages produce the same checking results by using both the old huge and the new huge.
