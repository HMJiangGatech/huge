#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
List SFGen(int d0, int d){
    IntegerMatrix G(d,d);
    int i,j;
    double x;
    int *size_a = (int*) Calloc(d, int);
    int tmp;
    int total;

    for(i=0;i<(d0-1);i++){
        G(i+1,i) = 1;
        G(i,i+1) = 1;
    }
    G(d0-1,0) = 1;
    G(0,d0-1) = 1;

    for(i=0;i<d0;i++){
        size_a[i] = 2;
    }
    for(i=d0;i<d;i++){
        size_a[i] = 0;
    }
    total = 2*d0;

    for(i=d0;i<d;i++){
        GetRNGstate();
        x = (double) total*unif_rand();
        PutRNGstate();
        tmp = 0;
        j = 0;
        while(tmp<x&&j<i){
            tmp += size_a[j];
            j++;
        }
        j = j-1;
        G(j,i) = 1;
        G(i,j) = 1;
        total+=2;
        size_a[j]++;
        size_a[i]++;
    }
    Free(size_a);

    return List::create(
      _["G"] = G
    );
}
