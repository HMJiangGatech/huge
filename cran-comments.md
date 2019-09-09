## Test environments
* Mac OS X, R 3.6.1
* ubuntu 16.04, R 3.6.1
* windows 10, R 3.6.1

## News
Compared to 1.3.2, we fixed a minor bug of using symmetric representation. The api and the others remains the same. We also fixed a major bug that makes the package very memory consuming. 

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:
  
* checking installed package size ... NOTE
  installed size is 9.8Mb
  sub-directories of 1Mb or more:
    libs 8.1Mb
  
  This is a standard NOTE when the C/C++ code is included.