
data {
  int<lower=1> Ndat; // number of datasets
  int<lower=1> N[Ndat]; //individuals in each dataset
  int<lower=1> rep_d[Ndat];//number of stool samples per individual
  int<lower=1> rep_s[Ndat];//number of slides per stool sample
  int<lower=1> N_kk;//total number of KK datapoints
  int<lower=-1> y_kk[N_kk];//vector of KK measurements (repeat slides, repeat day, individuals, datasets)
  int<lower=1> N_ind;//total number of individuals
  int<lower=1> N_ind_day;//total number of measurement days
  int<lower=1> N_cca_non;//sum of #individuals*#CCA interpretations
}


parameters {
  //parameters that are unrelated between settings

  real<lower=0> pop_var[Ndat];
  simplex[2] theta[Ndat]; 
  real<lower=0.2> mean_inf[Ndat];
  //day variation parameters
  real<lower=0> d_var_0;
  real<lower=0> d_var_1;
  //nb parameter, variation between repeated slides
  real<lower=0> count_var_0;
  // real count_var_loc_var;
  // vector[Ndat] count_var_loc;
  
  //individual level parameters
  vector<lower=0>[N_ind] lambda;
  vector[N_ind_day] log_d_var_ind;
}


transformed parameters{


  vector[2] alpha[Ndat];
  //general transformed parameters
  real<lower=0> non_inf_mean;
  real<lower=0> non_inf_var;
  real<lower=0> min_infect;
  real<lower=0> sigma[N_ind_day,2]; // scales of mixture components
  
  //you know
  real<lower=0> count_var[N_ind_day];
  //the ones that have to be initiated for each dataset
  // mixing proportions
  vector[N_ind_day] mu[2]; // locations of mixture components
  non_inf_var=1000;  
  non_inf_mean=0.001;
  min_infect=0.2;//5; does offset actually do anything

  mu[1,]= rep_vector(non_inf_mean,N_ind_day);
  
  {
    int pos1;
    int pos3;
    real temp;
    pos1=0;
    pos3=0;
    
    
    for(i in 1:Ndat){
    alpha[i,1]=1;
    alpha[i,2]=1;
      for(j in 1:rep_d[i]){
        for(l in 1:N[i]){
          temp=log(1+1/(lambda[pos1+l]+min_infect)+1/(d_var_1*(lambda[pos1+l]+min_infect)));
          mu[2,pos3+(j-1)*N[i]+l]=(lambda[pos1+l]+min_infect)*exp(log_d_var_ind[pos3+(j-1)*N[i]+l]*temp-temp^2/2);
          count_var[pos3+(j-1)*N[i]+l]=count_var_0*(lambda[pos1+l]+min_infect);//+(lambda[pos1+l]+min_infect)*count_var_1;
          sigma[pos3+(j-1)*N[i]+l,2]=count_var[pos3+(j-1)*N[i]+l];
          sigma[pos3+(j-1)*N[i]+l,1]=non_inf_var;
        }
      }
      pos3=pos3+rep_d[i]*N[i];
      pos1= pos1 + N[i];
    }
  }
  
}
model{
  real ps[2]; // temp for log component densities
  int pos2;
  int pos3;
  int indi2;
  pos2=1;
  pos3=0;
  indi2=0;
  pop_var~normal(0.5,0.5);
  mean_inf~normal(30,20);
  d_var_0~normal(0,3);
  d_var_1~normal(0,3);
  count_var_0~normal(1,2);

  

  log_d_var_ind~normal(0,1);
  
  for(i in 1:Ndat){

    segment(lambda,pos2,N[i])~gamma((mean_inf[i]-min_infect)*pop_var[i],pop_var[i]);

    
    for(n in 1:N[i]){
      
      for(k in 1:2){
        ps[k]=log(theta[i,k]);
        for(j in 1:rep_d[i]){

          for(l in 1:rep_s[i]){
            if(y_kk[indi2+(n-1)*rep_s[i]*rep_d[i]+(j-1)*rep_s[i]+l]==-1){}
            else{
              ps[k] = ps[k] +  neg_binomial_2_lpmf(y_kk[indi2+(n-1)*rep_s[i]*rep_d[i]+(j-1)*rep_s[i]+l] | mu[k,pos3+(j-1)*N[i]+n], sigma[pos3+(j-1)*N[i]+n,k]);
            }
          }
        }
      }
      target += log_sum_exp(ps);
    }
    indi2=indi2+rep_d[i]*rep_s[i]*N[i];
    pos3=pos3+rep_d[i]*N[i];
    pos2=pos2+N[i];
  }
}
generated quantities{
  
}
