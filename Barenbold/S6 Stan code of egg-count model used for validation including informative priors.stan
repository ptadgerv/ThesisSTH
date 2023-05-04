data {
  int<lower=0> N;
  int<lower=0> Ns;
  int<lower=1,upper=Ns> school[N];
  int<lower=0> epsA[N];
  int<lower=0> epsB[N];

}
transformed data{
  int<lower=0> Ndat=Ns;
}
parameters {

  real<lower=0> d_var_0; 
  real<lower=0> d_var_1;
  real<lower=0> count_var_1;
  simplex[2] theta[Ns]; 
  real<lower=0> pop_var[Ns];
  real<lower=0> mean_inf[Ns];
  
  //individual level parameters
  vector<lower=0>[N] lambda;
  vector[N] log_d_var_ind;
}
transformed parameters{

    //general transformed parameters
  real<lower=0> count_var[N];
  real<lower=0> non_inf_mean;
  real<lower=0> non_inf_var;
  real<lower=0> min_infect;
  real<lower=0> sigma[N,2]; // scales of mixture components
  vector<lower=0>[N] mu[2]; // locations of mixturecomponents
  real<lower=0> temp;
  non_inf_var=1000;  
  non_inf_mean=0.001;
  min_infect=0.0;//try without min_infect to avoid problems
  
  
 //changed to invert
  mu[1,]= rep_vector(non_inf_mean,N);
    
  for(l in 1:N){
    temp=log(d_var_0+1/(lambda[l]+min_infect)+1/(d_var_1*(lambda[l]+min_infect)));
    mu[2,l]=(lambda[l]+min_infect)*exp(log_d_var_ind[l]*temp-temp^2/2);
    count_var[l]=count_var_1;//mu[2,l]*count_var_1;
    sigma[l,1]=non_inf_var;
    sigma[l,2]=count_var[l];
    }
}
model{
  real ps[2]; // temp for log component densities

  pop_var~normal(0.5,0.5);
  mean_inf~normal(30,20);
  d_var_0~normal(0.75,0.05);
  d_var_1~normal(2.5,1);
  count_var_1~normal(2.5,0.1);

  log_d_var_ind~normal(0,1);
  for(n in 1:N){
      lambda[n]~gamma((mean_inf[school[n]]-min_infect)*pop_var[school[n]],pop_var[school[n]]); 
      for(k in 1:2){
        ps[k]=log(theta[school[n],k]);
        ps[k] = ps[k] +  neg_binomial_2_lpmf(epsA[n] | mu[k,n], sigma[n,k]);
        ps[k] = ps[k] +  neg_binomial_2_lpmf(epsB[n] | mu[k,n], sigma[n,k]); 
      }
      target += log_sum_exp(ps);
  }

}
