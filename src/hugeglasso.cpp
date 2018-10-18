#include "math.h"
#include <vector>
#include <omp.h>
#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;
using namespace Eigen;
//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::plugins(openmp)]

void hugeglasso_sub(Eigen::MatrixXd &S, Eigen::MatrixXd &W, Eigen::MatrixXd &T, int d, double ilambda, int &df, bool scr);

//[[Rcpp::export]]
List hugeglasso(Eigen::Map<Eigen::MatrixXd> S, NumericVector lambda, bool scr, bool verbose, bool cov_output)
{
    unsigned int nlambda = lambda.size();
    int nfeatures = S.rows();
    auto S_diag = S.diagonal().array();
    int &d = nfeatures;

    bool zero_sol = true;

    // define result list: loglik
    NumericVector loglik(nlambda,-nfeatures), sparsity(nlambda);
    IntegerVector df(nlambda);
    List result;
    result["loglik"] = loglik;
    result["sparsity"] = sparsity;
    result["df"] = df;

    std::vector<Eigen::MatrixXd *> tmp_icov_p, tmp_cov_p, tmp_path_p;
    for(int i = nlambda-1; i >= 0; i--)
    {
        // pre-screening by z
        vector<int> z;
        for (size_t row_i = 0; row_i < d; row_i++) {
            int break_flag = 0;
            for (size_t col_i = 0; col_i < d; col_i++) {
                if(break_flag > 1) break;
                if(S(row_i,col_i) > lambda[i] or S(row_i,col_i) < -lambda[i]) break_flag++;
            }
            if(break_flag > 1) z.push_back(row_i);
        }
        int q = z.size();
        MatrixXd sub_S(q,q), sub_W(q,q), sub_T(q,q);
        int sub_df=0;

        //#pragma omp parallel for
        for (size_t ii = 0; ii < q; ii++) {
            for (size_t jj = 0; jj < q; jj++) {
                sub_S(ii,jj) = S(z[ii],z[jj]);
                if(zero_sol) {
                    sub_W(ii,jj) = S(z[ii],z[jj]);
                    sub_T(ii,jj) = ii==jj ? 1 : 0;
                }
                else {
                    sub_W(ii,jj) = (*(tmp_cov_p.back()))(z[ii],z[jj]);
                    sub_T(ii,jj) = (*(tmp_icov_p.back()))(z[ii],z[jj]);
                }
            }
        }

        if(q>0)
        {
            if(verbose){
              if(scr)
                Rcout << "\rConducting the graphical lasso (glasso) wtih lossy screening....in progress: " << floor(100*(1-1.*i/nlambda))<<"%";
              if(!scr)
                Rcout << "\rConducting the graphical lasso (glasso) wtih lossless screening....in progress: " << floor(100*(1-1.*i/nlambda))<<"%";
            }

            hugeglasso_sub(sub_S, sub_W, sub_T, q, lambda[i], sub_df, scr);
            zero_sol = false;
        }
        if(q == 0) zero_sol = true;
        // update result list
        tmp_path_p.push_back(new Eigen::MatrixXd(d,d));
        tmp_icov_p.push_back(new Eigen::MatrixXd(d,d));
        tmp_cov_p.push_back(new Eigen::MatrixXd(d,d));

        Eigen::MatrixXd *tmp_icov, *tmp_cov, *tmp_path;
        tmp_icov=tmp_icov_p.back();
        tmp_cov=tmp_cov_p.back();
        tmp_path=tmp_path_p.back();
        tmp_icov->setZero();
        tmp_icov->diagonal() = 1/(S_diag+double(lambda[i]));
        tmp_cov->setZero();
        tmp_cov->diagonal() = S_diag + double(lambda[i]);
        tmp_path->setZero();

        if(!zero_sol)
        {
            //#pragma omp parallel for
            for (size_t ii = 0; ii < q; ii++) {
                for (size_t jj = 0; jj < q; jj++) {
                    (*tmp_icov)(z[ii],z[jj]) = sub_T(ii,jj);
                    (*tmp_cov)(z[ii],z[jj]) = sub_W(ii,jj);
                    (*tmp_path)(z[ii],z[jj]) = sub_T(ii,jj)==0 ? 0 : 1;
                    (*tmp_path)(z[ii],z[jj]) = ii==jj ? 0 : (*tmp_path)(z[ii],z[jj]);
                }
            }
            sparsity[i] = 1.0*sub_df/d/(d-1);
            df[i] = sub_df/2;
            loglik[i] = log(sub_T.determinant()) - (sub_T * sub_S).diagonal().sum() - (d-q);
        }
    }

    List path, icov, cov;
    for(int i = 0; i < nlambda; i++){
        path.push_back(*(tmp_path_p[nlambda-1-i]));
        icov.push_back(*(tmp_icov_p[nlambda-1-i]));
        if(cov_output) cov.push_back(*(tmp_cov_p[nlambda-1-i]));
    }
    result["path"] = path;
    result["icov"] = icov;
    if(cov_output) result["cov"] = cov;
    return result;
}

void hugeglasso_sub(Eigen::MatrixXd &S, Eigen::MatrixXd &W, Eigen::MatrixXd &T, int d, double ilambda, int &df, bool scr)
{
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
    //#pragma omp parallel for
    for(i=0;i<d;i++){

        W(i, i) = S(i, i) + ilambda; //The diagonal elements are set optimal
        size_a[i] = 0;
        tmp1 = T(i, i);
        T(i, i) = 0;

        for(j=0;j<d;j++){
            if(scr)
                if(fabs(S(j, i)) <= ilambda){
                  idx_i(j, i) = -1;
                  T(j, i) = 0;
                  continue;
                }

            if(T(j, i)!=0){
                idx_a(size_a[i], i) = j; //initialize the active set
                size_a[i]++;
                idx_i(j, i) = -1; //initialize the inactive set
                T(j, i) = -T(j, i)/tmp1;
            }
            else idx_i(j, i) = 1;
        }
        idx_i(i, i) = -1;
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
}
