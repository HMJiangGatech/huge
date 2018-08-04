#include "math.h"
#include <omp.h>
#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::plugins(openmp)]
//[[Rcpp::export]]
List hugeglasso(Eigen::MatrixXd &S, Eigen::MatrixXd &W, Eigen::MatrixXd &T, int d, double ilambda)
{
    int df = 0;

    int i,j,k; //initialize indices
    int rss_idx,w_idx;

    int gap_int;
    double gap_ext,gap_act;
    double thol_act = 1e-4;
    double thol_ext = 1e-4;

    int MAX_ITER_EXT = 100;
    int MAX_ITER_INT = 10000;
    int MAX_ITER_ACT = 10000;
    int iter_ext,iter_int,iter_act;



    Eigen::MatrixXi idx_a(d, d); // active set
    Eigen::MatrixXi idx_i(d, d); // The set possibly can join active set
    int *size_a = (int*) malloc(d*sizeof(int)); //sizes of active sets
    double *w1 = (double*) malloc(d*sizeof(double));
    double *ww = (double*) malloc(d*sizeof(double));


    int size_a_prev; //original size of the active set
    int junk_a; //the number of variables returning to the inactive set from the active set

    double r; //partial residual
    double tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;

    //Given the initial input W and T, recover inital solution for each individual lasso
    for(i=0;i<d;i++){

        W(i, i) = S(i, i) + ilambda; //The diagonal elements are set optimal
        size_a[i] = 0;
        tmp1 = T(i, i);
        T(i, i) = 0;
        idx_i(i, i) = -1;
        for(j=0;j<i;j++){
            if(T(j, i)!=0){
                idx_a(size_a[i], i) = j; //initialize the active set
                size_a[i]++;
                idx_i(j, i) = -1; //initialize the inactive set
                T(j, i) = -T(j, i)/tmp1;
            }
            else idx_i(j, i) = 1;
        }
        for(j=i+1;j<d;j++){
            if(T(j, i)!=0){
                idx_a(size_a[i], i) = j; //initialize the active set
                size_a[i]++;
                idx_i(j, i) = -1; //initialize the inactive set
                T(j, i) = -T(j, i)/tmp1;
            }
            else idx_i(j, i) = 1;
        }
    }

    gap_ext = 1;
    iter_ext = 0;
    while(gap_ext>thol_ext && iter_ext < MAX_ITER_EXT) //outer loop
    {
        tmp1 = 0;
        tmp6 = 0;
        tmp5 = 0;
        for(i=0;i<d;i++)
        {


            gap_int = 1;
            iter_int = 0;

            for(j=0;j<d;j++)
                ww[j] = T(j, i);

            while(gap_int!=0 && iter_int<MAX_ITER_INT)
            {
                size_a_prev = size_a[i];
                for(j=0;j<d;j++)
                {
                    if(idx_i(j, i)!=-1)
                    {

                        r = S(j, i);
                        for(k=0;k<size_a[i];k++)
                        {
                            rss_idx = idx_a(k, i);
                            r = r - W(rss_idx, j)*T(rss_idx, i);
                        }
                        if(r>ilambda)
                        {
                            w1[j] = (r - ilambda)/W(j, j);
                            idx_a(size_a[i], i) = j;
                            size_a[i] = size_a[i] + 1;
                            idx_i(j, i) = -1;
                        }

                        else if(r<-ilambda)
                        {
                            w1[j] = (r + ilambda)/W(j, j);
                            idx_a(size_a[i], i) = j;
                            size_a[i] = size_a[i] + 1;
                            idx_i(j, i) = -1;
                        }

                        else w1[j] = 0;

                        T(j, i) = w1[j];
                    }
                }
                gap_int = size_a[i] - size_a_prev;

                gap_act = 1;
                iter_act = 0;

                while(gap_act>thol_act && iter_act < MAX_ITER_ACT)
                {
                    tmp3 = 0;
                    tmp4 = 0;
                    for(j=0;j<size_a[i];j++)
                    {
                        w_idx = idx_a(j, i);
                        if(w_idx!=-1)
                        {
                            //tmp_a = w_idx*d;
                            r = S(w_idx, i) + T(w_idx, i)*W(w_idx, w_idx);
                            for(k=0;k<size_a[i];k++)
                            {
                                rss_idx = idx_a(k, i);
                                r = r - W(rss_idx, w_idx)*T(rss_idx, i);
                            }

                            if(r>ilambda){
                                w1[w_idx] = (r - ilambda)/W(w_idx, w_idx);
                                tmp4 += w1[w_idx];
                            }


                            else if(r<-ilambda){
                                w1[w_idx] = (r + ilambda)/W(w_idx, w_idx);
                                tmp4 -= w1[w_idx];
                            }

                            else w1[w_idx] = 0;

                            tmp3 = tmp3 + fabs(w1[w_idx] - T(w_idx, i));

                            T(w_idx, i) = w1[w_idx];
                        }
                    }
                    gap_act = tmp3/tmp4;
                    iter_act++;
                }

                //move the false active variables to the inactive set

                junk_a = 0;
                for(j=0;j<size_a[i];j++){
                    w_idx = idx_a(j, i);
                    if(w1[w_idx]==0){
                        junk_a++;
                        idx_i(w_idx, i) = 1;
                        idx_a(j, i) = -1;
                    }
                    else idx_a(j-junk_a, i) = w_idx;
                }
                size_a[i] = size_a[i] - junk_a;
                iter_int++;
            }

            //update W Beta
            Eigen::MatrixXd temp = W.transpose()*T.col(i);
            for(j=0;j<i;j++){
              W(j, i) = temp(j, 0);
              W(i, j) = temp(j, 0);
            }
            for(j=i+1;j<d;j++){
              W(j, i) = temp(j, 0);
              W(i, j) = temp(j, 0);
            }


            for(j=0;j<d;j++)
                tmp5 = tmp5 + fabs(ww[j]-T(j, i));
            tmp6 = tmp6 + tmp4;
        }
        gap_ext = tmp5/tmp6;
        //printf("%g\n",gap_ext);
        iter_ext++;
    }

    for(i=0;i<d;i++) //Compute the final T
    {
        tmp2 = 0;
        tmp2 = W.col(i).transpose()*T.col(i) - W(i, i)*T(i, i);

        tmp1 = 1/(W(i, i)-tmp2);
        T.col(i) *= -tmp1;
        T(i, i) = tmp1;
    }
    for(i=0;i<d;i++)
        df += size_a[i];

    free(size_a);
    free(w1);
    free(ww);
    return List::create(Named("W") = W, Named("T") = T, Named("df") = df);
}
