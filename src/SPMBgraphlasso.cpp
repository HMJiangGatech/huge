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
List SPMBgraphlasso(Eigen::Map<Eigen::MatrixXd> data, NumericVector lambda, int nlambda, int d, bool scr, IntegerMatrix idx_scr, int nscr)
{
    Eigen::ArrayXd r, grad, w1, Y, XX;
    Eigen::ArrayXXd X;
    int n = data.rows();
    int maxdf = 0;
    maxdf = (n < d ? n : d)*d;

    NumericVector x(d*maxdf*nlambda);
    IntegerVector col_cnz(d+1);
    IntegerVector row_idx(d*maxdf*nlambda);


    X = data;
    XX.resize(d);
    for (int j = 0; j < d; j++)
      XX[j] = (X.col(j)*X.col(j)).sum()/n;
    double prec = 1e-4;
    int max_iter = 1000;
    int cnz = 0;
    for(int m=0;m<d;m++)
    {
      grad.resize(d);
      grad.setZero();
      w1.resize(d);
      w1.setZero();
      r.resize(n);
      r.setZero();
      Y = X.col(m);

      std::vector<int> actset_indcat(d, 0);
      std::vector<int> actset_indcat_aux(d, 0);
      std::vector<int> actset_idx;
      actset_idx.clear();
      std::vector<double> gr(d, 0);

      r = Y;

      double v = 0.0;
      for (int i = 0; i < n; i++) {
        double pred = w1.matrix().dot(X.row(i).matrix());
        v += (Y[i] - pred) * (Y[i] - pred);
      }
      v = v / n;
      double dev_thr = fabs(v) * prec;
      for(int i = 0; i < d; i++) gr[i] = grad[i];

      double temp = 0;
      int flag1 = 0;
      int flag2 = 1;
      for(int i=0;i<nlambda;i++)
      {
        double ilambda = lambda[i];
        if(scr)
          for(int jj = 0; jj < nscr; jj++)
          {
            int j = idx_scr(m,jj);
            if (j==m) continue;
            if (actset_indcat[j] == 0)
            {
              temp = threshold(fabs(gr[j]), ilambda);
              if (fabs(temp) > 1e-8) actset_indcat[j] = 1;
            }
          }
        else
          for(int j = 0; j < d; j++)
          {
            if (j==m) continue;
            if (actset_indcat[j] == 0)
            {
              temp = threshold(fabs(gr[j]), ilambda);
              if (fabs(temp) > 1e-8) actset_indcat[j] = 1;
            }
          }

        int loopcnt_level_0 = 0;
        flag2 = 1;
        while (loopcnt_level_0 < max_iter)
        {
          loopcnt_level_0 += 1;

          // Step 1: First pass constructing active set
          bool new_active_idx = true;
          int loopcnt_level_1 = 0;
          if (flag1 * flag2 != 0)
          {
            loopcnt_level_1 = max_iter + 1;
            new_active_idx = true;
          }

          while (loopcnt_level_1 < max_iter)
          {
            loopcnt_level_1 += 1;
            bool terminate_loop_level_1 = false;

            for (int j = 0; j < m; j++)
            {
              if (actset_indcat[j] == 0) continue;

              double w1_old = w1[j];

              grad[j] = (r*X.col(j)).sum()/n;

              //w1_old = w1[j];
              double tmp = grad[j] + w1[j] * XX[j];
              w1[j] = threshold(tmp, ilambda) / XX[j];
              r = r - X.col(j) * (w1[j] - w1_old);
              double updated_coord = w1[j];

              if (updated_coord == w1_old) continue;

              if (actset_indcat_aux[j] == 0)
              {
                actset_idx.push_back(j);
                actset_indcat_aux[j] = 1;
              }
              tmp = w1_old - w1[j];
              double local_change = tmp * tmp * XX[j];
              if (local_change > dev_thr)
                terminate_loop_level_1 = true;
            }

            for (int j = m+1; j < d; j++)
            {
              if (actset_indcat[j] == 0) continue;

              double w1_old = w1[j];

              grad[j] = (r*X.col(j)).sum()/n;

              //w1_old = w1[j];
              double tmp = grad[j] + w1[j] * XX[j];
              w1[j] = threshold(tmp, ilambda) / XX[j];
              r = r - X.col(j) * (w1[j] - w1_old);
              double updated_coord = w1[j];

              if (updated_coord == w1_old) continue;

              if (actset_indcat_aux[j] == 0)
              {
                actset_idx.push_back(j);
                actset_indcat_aux[j] = 1;
              }
              tmp = w1_old - w1[j];
              double local_change = tmp * tmp * XX[j];
              if (local_change > dev_thr)
                terminate_loop_level_1 = true;
            }

            if (terminate_loop_level_1)
            {
              new_active_idx = true;
              break;
            }

            new_active_idx = false;
            for (int j = 0; j < m; j++)
              if (actset_indcat[j] == 0)
              {
                grad[j] = (r*X.col(j)).sum()/n;
                gr[j] = fabs(grad[j]);
                double temp = threshold(gr[j], ilambda);
                if (fabs(temp) > 1e-8) {
                  actset_indcat[j] = 1;
                  new_active_idx = true;
                }
              }

              for (int j = m+1; j < d; j++)
                if (actset_indcat[j] == 0)
                {
                  grad[j] = (r*X.col(j)).sum()/n;
                  gr[j] = fabs(grad[j]);
                  double temp = threshold(gr[j], ilambda);
                  if (fabs(temp) > 1e-8) {
                    actset_indcat[j] = 1;
                    new_active_idx = true;
                  }
                }

            if (!new_active_idx) break;
          }

          flag1 = 1;

          if (!new_active_idx) break;

          // Step 2 : active set minimization
          // on the active coordinates
          loopcnt_level_1 = 0;
          while (loopcnt_level_1 < max_iter)
          {
            loopcnt_level_1 += 1;

            bool terminate_loop_level_1 = true;
            for (unsigned int j = 0; j < actset_idx.size(); j++) {
              int idx = actset_idx[j];
              double w1_old = w1[idx];

              // compute thresholded coordinate
              grad[idx] = (r*X.col(idx)).sum()/n;
              w1_old = w1[idx];
              double tmp = grad[idx] + w1[idx] * XX[idx];
              w1[idx] = threshold(tmp, ilambda) / XX[idx];



              r = r - X.col(idx) * (w1[idx] - w1_old);
              double updated_coord = w1[idx];

              if (w1_old == updated_coord) continue;

              tmp = w1_old - w1[idx];
              double local_change = tmp * tmp * XX[idx];
              if (local_change > dev_thr)
                terminate_loop_level_1 = false;
              // update gradient for idx
            }

            if (terminate_loop_level_1) {
              flag2 = 0;
              break;
            }
          }

        }

        for(unsigned int j=0;j<actset_idx.size();j++)
        {
          int w_idx = actset_idx[j];
          x[cnz] = w1[w_idx];
          row_idx[cnz] = i*d+w_idx;
          cnz++;
        }
      }
      col_cnz[m+1]=cnz;
    }
    return List::create(
      _["col_cnz"] = col_cnz,
      _["row_idx"] = row_idx,
      _["x"] = x
    );
}


double threshold(double x, double thr) {
  if (x > thr)
    return x - thr;
  else if (x < -thr)
    return x + thr;
  else
    return 0;
}
