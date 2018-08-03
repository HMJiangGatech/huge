#include "R.h"
#include "math.h"

void RIC(double *X, int *dd, int *nn, int *r, int *nr, double *lambda_opt)
{
	int d,n;
	int i,j,k,m;
    int t;
    int tmp_r,tmp_i,tmp_j;
	
	d = dd[0];
	n = nn[0];
    t = nr[0];
    
	double lambda_min,lambda_max,tmp;
	
	lambda_min = 99999999;
	
    for(i=0;i<t;i++)
    {
        tmp_r = r[i];
        lambda_max = 0;
        for(j=0;j<d;j++)
            for(k=0;k<j;k++)
            {
                tmp = 0;
                for(m=0;m<(n-tmp_r);m++)
                    tmp = tmp + X[j*n+m+tmp_r]*X[k*n+m];
                for(m=(n-tmp_r);m<n;m++)
                    tmp = tmp + X[j*n+m-(n-tmp_r)]*X[k*n+m];
                tmp = fabs(tmp);
                if(tmp>lambda_max)
                    lambda_max = tmp;
            }
            for(k=j+1;k<d;k++)
            {
                tmp = 0;
                for(m=0;m<(n-tmp_r);m++)
                    tmp = tmp + X[j*n+m+tmp_r]*X[k*n+m];
                for(m=(n-tmp_r);m<n;m++)
                    tmp = tmp + X[j*n+m-(n-tmp_r)]*X[k*n+m];
                tmp = fabs(tmp);
                if(tmp>lambda_max)
                    lambda_max = tmp;
            }
        if(lambda_max<lambda_min)
            lambda_min = lambda_max;
    }
    lambda_opt[0] = lambda_min;
}
