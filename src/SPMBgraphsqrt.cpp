#include "math.h"
#include <vector>
#include <RcppEigen.h>
#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;
using namespace Eigen;
//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::plugins(openmp)]


double thresholdl1(double x, double thr);

//[[Rcpp::export]]
List SPMBgraphsqrt(Eigen::Map<Eigen::MatrixXd> data, NumericVector lambda, int nlambda, int d)
{

    Eigen::ArrayXd Xb, r, grad, w1, Y, XX, gr;
    Eigen::ArrayXXd X;
    Eigen::MatrixXd tmp_icov;
    tmp_icov.resize(d, d);
    tmp_icov.setZero();
    std::vector<Eigen::MatrixXd > tmp_icov_p;
    tmp_icov_p.clear();
    for(int i = 0; i < nlambda; i++)
      tmp_icov_p.push_back(tmp_icov);
    int n = data.rows();
    int maxdf = 0;
    maxdf = (n < d ? n : d);
    NumericVector x(d*maxdf*nlambda);
    IntegerVector col_cnz(d+1);
    IntegerVector row_idx(d*maxdf*nlambda);
    X = data;
    XX.resize(d);
    for (int j = 0; j < d; j++)
      XX[j] = (X.col(j)*X.col(j)).sum()/n;
    double prec = 1e-4;
    int max_iter = 1000;
    int num_relaxation_round = 3;
	  int cnz = 0;
    for(int m=0;m<d;m++)
    {
      Xb.resize(n);
      Xb.setZero();
      grad.resize(d);
      grad.setZero();
      gr.resize(d);
      gr.setZero();
      w1.resize(d);
      w1.setZero();
      r.resize(n);
      r.setZero();
      Y = X.col(m);

      Eigen::ArrayXd Xb_master(n);
      Eigen::ArrayXd w1_master(n);
      std::vector<int> actset_indcat(d, 0);
      std::vector<int> actset_indcat_master(d, 0);
      std::vector<int> actset_idx;
      std::vector<double> old_coef(d, 0);
      std::vector<double> grad(d, 0);
      std::vector<double> grad_master(d, 0);

      double a = 0, g = 0, L = 0, sum_r2 = 0;
      double tmp_change = 0, local_change = 0;

      r = Y - Xb;
      sum_r2 = r.matrix().dot(r.matrix());
      L = sqrt(sum_r2 / n);

      double dev_thr = fabs(L) * prec;

      //cout<<dev_thr<<endl;


      for(int i = 0; i < d; i++)
      {
        grad[i] = (r * X.col(i)).sum() / (n*L);
      }
      for(int i = 0; i < d; i++) gr[i] = abs(grad[i]);
      w1_master = w1;
      Xb_master = Xb;
      for (int i = 0; i < d; i++) grad_master[i] = gr[i];

      std::vector<double> stage_lambdas(d, 0);

      for(int i=0;i<nlambda;i++)
      {
        w1 = w1_master;
        Xb = Xb_master;

        for (int j = 0; j < d; j++)
        {
          gr[j] = grad_master[j];
          actset_indcat[j] = actset_indcat_master[j];
        }

        // init the active set
        double threshold;
        if (i > 0)
          threshold = 2 * lambda[i] - lambda[i - 1];
        else
          threshold = 2 * lambda[i];

        for (int j = 0; j < m; ++j)
        {
          stage_lambdas[j] = lambda[i];

          if (gr[j] > threshold) actset_indcat[j] = 1;
        }
        for (int j = m+1; j < d; ++j)
        {
          stage_lambdas[j] = lambda[i];

          if (gr[j] > threshold) actset_indcat[j] = 1;
        }
        stage_lambdas[m] = lambda[i];
        r = Y - Xb;
        sum_r2 = r.matrix().dot(r.matrix());
        L = sqrt(sum_r2 / n);
        // loop level 0: multistage convex relaxation
        int loopcnt_level_0 = 0;
        int idx;
        double old_w1, updated_coord;
        while (loopcnt_level_0 < num_relaxation_round)
        {
          loopcnt_level_0++;

          // loop level 1: active set update
          int loopcnt_level_1 = 0;
          bool terminate_loop_level_1 = true;
          while (loopcnt_level_1 < max_iter)
          {
            loopcnt_level_1++;
            terminate_loop_level_1 = true;

            for (int j = 0; j < d; j++) old_coef[j] = w1[j];

            // initialize actset_idx
            actset_idx.clear();
            for (int j = 0; j < m; j++)
              if (actset_indcat[j])
              {
                g = 0.0;
                a = 0.0;

                double tmp;

                sum_r2 = r.matrix().dot(r.matrix());
                L = sqrt(sum_r2 / n);

                Eigen::ArrayXd wXX  = (1 - r*r/sum_r2) * X.col(j) * X.col(j);
                g = (wXX * w1[j] + r * X.col(j)).sum()/(n*L);
                a = wXX.sum()/(n*L);

                tmp = w1[j];
                w1[j] = thresholdl1(g, stage_lambdas[j]) / a;

                tmp = w1[j] - tmp;
                // Xb += delta*X[idx*n]
                Xb = Xb + tmp * X.col(j);

                sum_r2 = 0.0;
                // r -= delta*X
                r = r - tmp * X.col(j);

                sum_r2 = r.matrix().dot(r.matrix());
                L = sqrt(sum_r2 / n);

                updated_coord = w1[j];

                if (fabs(updated_coord) > 0) actset_idx.push_back(j);
              }

            for (int j = m+1; j < d; j++)
              if (actset_indcat[j])
              {
                  g = 0.0;
                  a = 0.0;

                  double tmp;

                  sum_r2 = r.matrix().dot(r.matrix());
                  L = sqrt(sum_r2 / n);

                  Eigen::ArrayXd wXX  = (1 - r*r/sum_r2) * X.col(j) * X.col(j);
                  g = (wXX * w1[j] + r * X.col(j)).sum()/(n*L);
                  a = wXX.sum()/(n*L);

                  tmp = w1[j];
                  w1[j] = thresholdl1(g, stage_lambdas[j]) / a;

                  tmp = w1[j] - tmp;
                  // Xb += delta*X[idx*n]
                  Xb = Xb + tmp * X.col(j);

                  sum_r2 = 0.0;
                  // r -= delta*X
                  r = r - tmp * X.col(j);

                  sum_r2 = r.matrix().dot(r.matrix());
                  L = sqrt(sum_r2 / n);

                  updated_coord = w1[j];

                  if (fabs(updated_coord) > 0) actset_idx.push_back(j);
                }

              // loop level 2: proximal newton on active set
            int loopcnt_level_2 = 0;
            bool terminate_loop_level_2 = true;
            while (loopcnt_level_2 < max_iter)
            {
              loopcnt_level_2++;
              terminate_loop_level_2 = true;

              for (unsigned int k = 0; k < actset_idx.size(); k++)
              {
                  idx = actset_idx[k];

                  old_w1 = w1[idx];
                  g = 0.0;
                  a = 0.0;

                  double tmp;

                  sum_r2 = r.matrix().dot(r.matrix());
                  L = sqrt(sum_r2 / n);

                  Eigen::ArrayXd wXX  = (1 - r*r/sum_r2) * X.col(idx) * X.col(idx);
                  g = (wXX * w1[idx] + r * X.col(idx)).sum()/(n*L);
                  a = wXX.sum()/(n*L);

                  tmp = w1[idx];
                  w1[idx] = thresholdl1(g, stage_lambdas[idx]) / a;

                  tmp = w1[idx] - tmp;
                  // Xb += delta*X[idx*n]
                  Xb = Xb + tmp * X.col(idx);

                  sum_r2 = 0.0;
                  // r -= delta*X
                  r = r - tmp * X.col(idx);

                  sum_r2 = r.matrix().dot(r.matrix());
                  L = sqrt(sum_r2 / n);

                  updated_coord = w1[idx];
                  tmp_change = old_w1 - w1[idx];
                  double a =  (X.col(idx) * X.col(idx) * (1 - r * r/(L*L*n))).sum()/(n*L);
                  local_change = a * tmp_change * tmp_change / (2 * L * n);
                  if (local_change > dev_thr)
                    terminate_loop_level_2 = false;
                }
              if (terminate_loop_level_2)
                break;
            }

            terminate_loop_level_1 = true;
              // check stopping criterion 1: fvalue change
            for (unsigned int k = 0; k < actset_idx.size(); ++k)
            {
              idx = actset_idx[k];
              tmp_change = old_w1 - w1[idx];
              double a =  (X.col(idx) * X.col(idx) * (1 - r * r/(L*L*n))).sum()/(n*L);
              local_change = a * tmp_change * tmp_change / (2 * L * n);
              if (local_change > dev_thr)
                terminate_loop_level_1 = false;
            }

            r = Y - Xb;
            sum_r2 = r.matrix().dot(r.matrix());
            L = sqrt(sum_r2 / n);

            if (terminate_loop_level_1)
              break;


              // check stopping criterion 2: active set change
            bool new_active_idx = false;
            for (int k = 0; k < m; k++)
              if (actset_indcat[k] == 0)
              {
                grad[idx] = (r * X.col(idx)).sum() / (n*L);
                //cout<<grad[idx];
                gr[k] = fabs(grad[k]);
                if (gr[k] > stage_lambdas[k])
                {
                  actset_indcat[k] = 1;
                  new_active_idx = true;
                }
              }
            for (int k = m+1; k < d; k++)
              if (actset_indcat[k] == 0)
              {
                grad[idx] = (r * X.col(idx)).sum() / (n*L);
                //cout<<grad[idx]
                gr[k] = fabs(grad[k]);
                if (gr[k] > stage_lambdas[k])
                {
                  actset_indcat[k] = 1;
                  new_active_idx = true;
                }
              }
            if(!new_active_idx)
              break;
          }
          if (loopcnt_level_0 == 1)
          {
            for (int j = 0; j < d; j++)
            {
              w1_master[j] = w1[j];

              grad_master[j] = gr[j];
              actset_indcat_master[j] = actset_indcat[j];
            }

            for (int j = 0; j < n; j++) Xb_master[j] = Xb[j];
          }
        }
        for(unsigned int j=0;j<actset_idx.size();j++)
        {
          int w_idx = actset_idx[j];
          x[cnz] = w1[w_idx];
          row_idx[cnz] = i*d+w_idx;
          cnz++;
          //cout<<cnz<<"    ";
        }
        double tal = 0;
        Eigen::MatrixXd temp;
        temp.resize(n, 1);
        for(int j = 0; j < n; j++)
        {
          temp(j, 0) = 0;
          for(int k = 0; k < d; k++)
            temp(j, 0) += X.matrix()(j, k)*w1[k];
          temp(j, 0) = Y[j] - temp(j, 0);
        }
        //temp = Y.matrix() - X.matrix().transpose()*w1.matrix();
        for(int j = 0; j < n; j++)
          tal += temp(j, 0)*temp(j, 0);
        tal = sqrt(tal)/sqrt(n);

        tmp_icov = tmp_icov_p[i];
        tmp_icov(m, m) = pow(tal, -2);
        for(int j = 0; j < m; j++)
          tmp_icov(j, m) = -tmp_icov(m, m)*w1[j];
        for(int j = m+1; j < d; j++)
          tmp_icov(j, m) = -tmp_icov(m, m)*w1[j];
        tmp_icov_p[i] = tmp_icov;
      }
      col_cnz[m+1]=cnz;
    }
    for(int i = 0; i < nlambda; i++)
      tmp_icov_p[i] = (tmp_icov_p[i].transpose()+tmp_icov_p[i])/2;
	  return List::create(
	    _["col_cnz"] = col_cnz,
	    _["row_idx"] = row_idx,
	    _["x"] = x,
	    _["icov"] = tmp_icov_p
	  );
}


double thresholdl1(double x, double thr) {
  if (x > thr)
    return x - thr;
  else if (x < -thr)
    return x + thr;
  else
    return 0;
}
