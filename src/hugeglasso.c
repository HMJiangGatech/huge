#include "R.h"
#include "math.h"

void hugeglasso(double *S, double *W, double *T, int *dd, double *lambda, int *df) 
{
    
    double ilambda;
    int d;
    int d2;
    ilambda = lambda[0];
    d = dd[0];
    d2 = d*d;
    
    int i,j,k; //initialize indices
    int rss_idx,w_idx;
    int tmp_i;
    int tmp_j,tmp_a;
    
    
    int gap_int;
    double gap_ext,gap_act;
    double thol_act = 1e-4;
    double thol_ext = 1e-4;
    
    int MAX_ITER_EXT = 100;
    int MAX_ITER_INT = 10000;
    int MAX_ITER_ACT = 10000;
    int iter_ext,iter_int,iter_act;
    
    
    int *idx_a = (int*) malloc((d2)*sizeof(int)); //active sets
    int *idx_i = (int*) malloc((d2)*sizeof(int)); //inactive sets
    int *size_a = (int*) malloc(d*sizeof(int)); //sizes of active sets
    double *w1 = (double*) malloc(d*sizeof(double));
    double *ww = (double*) malloc(d*sizeof(double));
    
    int size_a_prev; //original size of the active set
    int junk_a; //the number of variables returning to the inactive set from the active set
    
    double r; //partial residual
    double tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
        
    //Given the initial input W and T, recover inital solution for each individual lasso
    for(i=0;i<d;i++){
        tmp_i = i*d;    
        W[tmp_i+i] = S[tmp_i+i] + ilambda; //The diagonal elements are set optimal
        size_a[i] = 0;
        tmp1 = T[tmp_i+i];
        T[tmp_i+i] = 0;
        idx_i[tmp_i+i] = -1;
        for(j=0;j<i;j++){
            if(T[tmp_i+j]!=0){
                idx_a[tmp_i+size_a[i]] = j; //initialize the active set
                size_a[i]++;
                idx_i[tmp_i+j] = -1; //initialize the inactive set
                T[tmp_i+j] = -T[tmp_i+j]/tmp1;
            }
            else idx_i[tmp_i+j] = 1;
        }
        for(j=i+1;j<d;j++){
            if(T[tmp_i+j]!=0){
                idx_a[tmp_i+size_a[i]] = j; //initialize the active set
                size_a[i]++;
                idx_i[tmp_i+j] = -1; //initialize the inactive set
                T[tmp_i+j] = -T[tmp_i+j]/tmp1;
            }
            else idx_i[tmp_i+j] = 1;
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
            
            tmp_i = i*d;
            gap_int = 1;
            iter_int = 0;
            
            for(j=0;j<d;j++)
                ww[j] = T[tmp_i+j];
            while(gap_int!=0 && iter_int<MAX_ITER_INT)
            { 
                size_a_prev = size_a[i];
                for(j=0;j<d;j++)
                {
                    if(idx_i[tmp_i+j]!=-1)
                    {
                        tmp_j = j*d;
                        r = S[tmp_i+j];
                        for(k=0;k<size_a[i];k++)
                        {
                            rss_idx = idx_a[tmp_i+k];
                            r = r - W[tmp_j+rss_idx]*T[tmp_i+rss_idx];
                        }
                        if(r>ilambda)
                        {
                            w1[j] = (r - ilambda)/W[tmp_j+j];
                            idx_a[tmp_i+size_a[i]] = j;
                            size_a[i] = size_a[i] + 1;
                            idx_i[tmp_i+j] = -1;
                            
                        }
                        
                        else if(r<-ilambda)
                        {
                            w1[j] = (r + ilambda)/W[tmp_j+j];
                            idx_a[tmp_i+size_a[i]] = j;
                            size_a[i] = size_a[i] + 1;
                            idx_i[tmp_i+j] = -1;
                        }
                        
                        else w1[j] = 0;
                        
                        T[tmp_i+j] = w1[j];
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
                        w_idx = idx_a[tmp_i+j];
                        if(w_idx!=-1)
                        {
                            tmp_a = w_idx*d;
                            r = S[tmp_i+w_idx] + T[tmp_i+w_idx]*W[tmp_a+w_idx];
                            for(k=0;k<size_a[i];k++)
                            {
                                rss_idx = idx_a[tmp_i+k];
                                r = r - W[tmp_a+rss_idx]*T[tmp_i+rss_idx];
                            }
                            
                            if(r>ilambda){
                                w1[w_idx] = (r - ilambda)/W[tmp_a+w_idx];
                                tmp4 += w1[w_idx];
                            }
                            
                            
                            else if(r<-ilambda){
                                w1[w_idx] = (r + ilambda)/W[tmp_a+w_idx];
                                tmp4 -= w1[w_idx];
                            }
                            
                            else w1[w_idx] = 0;
                            
                            tmp3 = tmp3 + fabs(w1[w_idx] - T[tmp_i+w_idx]);
                            
                            T[tmp_i+w_idx] = w1[w_idx];
                        }
                    }
                    gap_act = tmp3/tmp4;
                    iter_act++;
                }
                
                //move the false active variables to the inactive set
                
                junk_a = 0;
                for(j=0;j<size_a[i];j++){
                    w_idx = idx_a[tmp_i+j];
                    if(w1[w_idx]==0){
                        junk_a++;
                        idx_i[tmp_i+w_idx] = 1;
                        idx_a[tmp_i+j] = -1;
                    }
                    else idx_a[tmp_i+j-junk_a] = w_idx;
                }
                size_a[i] = size_a[i] - junk_a;
                iter_int++;
            }
            
            for(j=0;j<i;j++) //update W Beta
            {
                tmp2 = 0;
                tmp_j = j*d;
                for(k=0;k<i;k++)
                    tmp2 = tmp2 + T[tmp_i+k]*W[tmp_j+k];
                for(k=i+1;k<d;k++)
                    tmp2 = tmp2 + T[tmp_i+k]*W[tmp_j+k];
                W[tmp_i+j] = tmp2; 
                W[tmp_j+i] = tmp2;
                
            }
                        
            for(j=i+1;j<d;j++){
                tmp2 = 0;
                tmp_j = j*d;
                for(k=0;k<i;k++)
                    tmp2 = tmp2 + T[tmp_i+k]*W[tmp_j+k];
                for(k=i+1;k<d;k++)
                    tmp2 = tmp2 + T[tmp_i+k]*W[tmp_j+k];
                W[tmp_i+j] = tmp2; 
                W[tmp_j+i] = tmp2;
            }
            for(j=0;j<d;j++)
                tmp5 = tmp5 + fabs(ww[j]-T[tmp_i+j]);
            tmp6 = tmp6 + tmp4;
        }
        gap_ext = tmp5/tmp6;
        //printf("%g\n",gap_ext);
        iter_ext++;
    }
    for(i=0;i<d;i++) //Compute the final T
    {
        tmp2 = 0;
        for(j=0;j<i;j++)
            tmp2 = tmp2 + W[i*d+j]*T[i*d+j];
        for(j=i+1;j<d;j++)
            tmp2 = tmp2 + W[i*d+j]*T[i*d+j];
        
        tmp1 = 1/(W[i*d+i]-tmp2);
        T[i*d+i] = tmp1;
        for(j=0;j<i;j++)
            T[i*d+j] = -tmp1*T[i*d+j]; 
        for(j=i+1;j<d;j++)
            T[i*d+j] = -tmp1*T[i*d+j];
    } 
    for(i=0;i<d;i++)
        df[0] += size_a[i];
     
    free(idx_a);
    free(idx_i);
    free(size_a);
    free(w1);
    free(ww);
}
