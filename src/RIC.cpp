#include "math.h"
#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
double RIC(NumericMatrix &X, int d, int n, NumericVector &r, int t)
{
    int i,j,k,m;
    int tmp_r;

    double lambda_min,lambda_max,tmp;

    lambda_min = 99999999;

    for(i=0;i<t;i++)
    {
        tmp_r = r[i];
        lambda_max = 0;
        for(j=0;j<d;j++) {
            for(k=0;k<j;k++)
            {
                tmp = 0;
                for(m=0;m<(n-tmp_r);m++)
                    tmp = tmp + X(m+tmp_r, j)*X(m,k);
                for(m=(n-tmp_r);m<n;m++)
                    tmp = tmp + X(m-(n-tmp_r), j)*X(m,k);
                tmp = fabs(tmp);
                if(tmp>lambda_max)
                    lambda_max = tmp;
            }
            for(k=j+1;k<d;k++)
            {
                tmp = 0;
                for(m=0;m<(n-tmp_r);m++)
                    tmp = tmp + X(m+tmp_r, j)*X(m,k);
                for(m=(n-tmp_r);m<n;m++)
                    tmp = tmp + X(m-(n-tmp_r), j)*X(m,k);
                tmp = fabs(tmp);
                if(tmp>lambda_max)
                    lambda_max = tmp;
            }
        }
        if(lambda_max<lambda_min)
            lambda_min = lambda_max;
    }
    return lambda_min;
}
