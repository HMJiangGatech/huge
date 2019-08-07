## Test environments
* Mac OS X, R 3.6.1
* ubuntu 16.04, R 3.6.1
* windows 10, R 3.6.1

## News
Compared to 1.3.2, we fixed a minor bug of using symmetric representation. The api and the others remains the same. 

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 2 NOTE:

* checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
  'default5.eps'

  I have no idea where this file comes from, there is not such a file in huge package. Maybe it comes from the pdf of vignettes.  
  
* checking installed package size ... NOTE
  installed size is 9.8Mb
  sub-directories of 1Mb or more:
    libs 8.1Mb
  
  This is a standard NOTE when the C/C++ code is included.