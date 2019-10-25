#include "math.h"
#include <RcppEigen.h>
#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;
using namespace Eigen;
//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::plugins(openmp)]


double threshold(double x, double thr);

//[[Rcpp::export]]
List SPMBscr(Eigen::Map<Eigen::MatrixXd> S, NumericVector lambda, int nlambda, int d, int maxdf, IntegerMatrix idx_scr, int nscr)
{
    NumericVector x(d*maxdf*nlambda);
    IntegerVector col_cnz(d+1);
    col_cnz[0] = 0;
    IntegerVector row_idx(d*maxdf*nlambda);int m,i,j,k;
    int cnz;
    int w_idx,rss_idx,size_a,size_a_prev;

    double ilambda;
    double r,tmp1,tmp2;

    int iter_ext,iter_int;
    int gap_ext;
    double gap_int;

    double thol = 1e-4;
    int MAX_ITER = 10000;

    double *w0 = (double*) malloc(d*sizeof(double));
    double *w1 = (double*) malloc(d*sizeof(double));
    int *idx_a = (int*) malloc(nscr*sizeof(int)); //sizes of active sets
    int *idx_i = (int*) malloc(nscr*sizeof(int)); //sizes of active sets

    cnz = 0;

    for(m=0;m<d;m++)
    {
        size_a = 0;

        for(j=0;j<nscr;j++)
            idx_i[j] = idx_scr[m*nscr+j];

        for(j=0;j<d;j++)
            w0[j] = 0;

        for(i=0;i<nlambda;i++)
        {
            ilambda = lambda[i];
            gap_ext = 1;
            iter_ext = 0;
            while(iter_ext<MAX_ITER && gap_ext>0)
            {
                size_a_prev = size_a;
                for(j=0;j<nscr;j++)
                {
                    w_idx = idx_i[j];
                    if(w_idx!=-1)
                    {
                        r = S(m,w_idx);

                        for(k=0;k<size_a;k++)
                        {
                            rss_idx = idx_a[k];
                            r = r - S(w_idx,rss_idx)*w0[rss_idx];
                        }

                        if(r > ilambda)
                        {
                            w1[w_idx] = r - ilambda;
                            idx_a[size_a] = w_idx;
                            size_a++;
                            idx_i[j] = -1;
                        }

                        if(r < -ilambda)
                        {
                            w1[w_idx] = r + ilambda;
                            idx_a[size_a] = w_idx;
                            size_a++;
                            idx_i[j] = -1;
                        }
                        else w1[w_idx] = 0;
                        w0[w_idx] = w1[w_idx];
                    }
                }

                gap_ext = size_a - size_a_prev;

                gap_int = 1;
                iter_int = 0;
                while(gap_int>thol && iter_int<MAX_ITER)
                {
                    tmp1 = 0;
                    tmp2 = 0;
                    for(j=0;j<size_a;j++)
                    {
                        w_idx = idx_a[j];
                        r = S(m,w_idx) + w0[w_idx];

                        for(k=0;k<size_a;k++)
                        {
                            rss_idx = idx_a[k];
                            r = r - S(w_idx,rss_idx)*w0[rss_idx];
                        }

                        if(r > ilambda)
                        {
                            w1[w_idx] = r - ilambda;
                            tmp2 += fabs(w1[w_idx]);
                        }

                        else if(r <-ilambda)
                        {
                            w1[w_idx] = r + ilambda;
                            tmp2 += fabs(w1[w_idx]);
                        }

                        else w1[w_idx] = 0;

                        tmp1 += fabs(w1[w_idx] - w0[w_idx]);
                        w0[w_idx] = w1[w_idx];
                    }
                    gap_int = tmp1/tmp2;
                    iter_int++;
                }
                iter_ext++;
            }
            for(j=0;j<size_a;j++)
            {
                w_idx = idx_a[j];
                x[cnz] = w1[w_idx];
                row_idx[cnz] = i*d+w_idx;
                cnz++;
            }
        }
        col_cnz[m+1]=cnz;
    }
    free(w0);
    free(w1);
    free(idx_a);
    free(idx_i);
    return List::create(
      _["col_cnz"] = col_cnz,
      _["row_idx"] = row_idx,
      _["x"] = x
    );
}

//[[Rcpp::export]]
List SPMBgraph(Eigen::Map<Eigen::MatrixXd> S, NumericVector lambda, int nlambda, int d, int maxdf)
{
    NumericVector x(d*maxdf*nlambda);
    IntegerVector col_cnz(d+1);
    col_cnz[0] = 0;
    IntegerVector row_idx(d*maxdf*nlambda);
    int m,i,j,k;

    int cnz,junk_a;
    int w_idx,rss_idx,size_a;
    int size_a_prev;

    double ilambda;
    double tmp1,tmp2;
    double r;

    int iter_ext,iter_int;
    int gap_ext;
    double gap_int;

    double thol = 1e-4;
    int MAX_ITER = 10000;

    double *w0 = (double*) malloc(d*sizeof(double));
    double *w1 = (double*) malloc(d*sizeof(double));
    int *idx_a = (int*) malloc(d*sizeof(int)); //sizes of active sets
    int *idx_i = (int*) malloc(d*sizeof(int)); //sizes of active sets

    cnz = 0;

    for(m=0;m<d;m++)
    {
        idx_i[m] = 0;
        for(j=0;j<m;j++)
            idx_i[j] = 1;
        for(j=m+1;j<d;j++)
            idx_i[j] = 1;

        size_a = 0;
        for(j=0;j<d;j++)
            w0[j] = 0;

        for(i=0;i<nlambda;i++)
        {
            ilambda = lambda[i];
            gap_ext = 1;
            iter_ext = 0;
            while(gap_ext !=0 && iter_ext<MAX_ITER)
            {
                size_a_prev = size_a;
                for(j=0;j<d;j++)
                {
                    if(idx_i[j]==1)
                    {
                        r = S(m,j);
                        for(k=0;k<size_a;k++)
                        {
                            rss_idx = idx_a[k];
                            r = r - S(j,rss_idx)*w0[rss_idx];
                        }

                        if(r > ilambda)
                        {
                            w1[j] = r - ilambda;
                            idx_a[size_a] = j;
                            size_a++;
                            idx_i[j] = 0;
                        }
                        else if(r <-ilambda)
                        {
                            w1[j] = r + ilambda;
                            idx_a[size_a] = j;
                            size_a++;
                            idx_i[j] = 0;
                        }

                        else w1[j] = 0;

                        w0[j] = w1[j];
                    }
                }

            gap_ext = size_a - size_a_prev;

            gap_int = 1;
            iter_int = 0;
            while(gap_int>thol && iter_int<MAX_ITER)
            {
                tmp1 = 0;
                tmp2 = 0;
                for(j=0;j<size_a;j++)
                {
                    //if(w_idx!=-1)
                    {
                        w_idx = idx_a[j];
                        r = S(m,w_idx) + w0[w_idx];

                        for(k=0;k<size_a;k++)
                        {
                            rss_idx = idx_a[k];
                            r = r - S(w_idx,rss_idx)*w0[rss_idx];
                        }

                        if(r > ilambda)
                        {
                            w1[w_idx] = r - ilambda;
                            tmp2 += fabs(w1[w_idx]);
                        }
                        else if(r <-ilambda){
                            w1[w_idx] = r + ilambda;
                            tmp2 += fabs(w1[w_idx]);
                        }

                        else w1[w_idx] = 0;

                        tmp1 += fabs(w1[w_idx] - w0[w_idx]);
                        w0[w_idx] = w1[w_idx];
                    }
                }
                gap_int = tmp1/tmp2;
                iter_int++;
            }
            junk_a = 0;
            for(j=0;j<size_a;j++)
            {
                w_idx = idx_a[j];
                if(w1[w_idx]==0){
                    junk_a++;
                    idx_i[w_idx] = 1;
                    //idx_a[j] = -1;
                }
                else idx_a[j-junk_a] = w_idx;
            }
            size_a = size_a - junk_a;
            iter_ext++;
            }
            for(j=0;j<size_a;j++)
            {
                w_idx = idx_a[j];
                x[cnz] = w1[w_idx];
                row_idx[cnz] = i*d+w_idx;
                cnz++;
            }
        }
        col_cnz[m+1]=cnz;
    }
    free(w0);
    free(w1);
    free(idx_a);
    free(idx_i);
    return List::create(
    _["col_cnz"] = col_cnz,
    _["row_idx"] = row_idx,
    _["x"] = x
    );
}
