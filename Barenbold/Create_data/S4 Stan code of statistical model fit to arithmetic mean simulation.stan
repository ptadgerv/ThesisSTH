data{
 int<lower=1> N_ind; //number of datapoints
 vector[N_ind] inf_eps;//epg infection intensity of each dataset
 vector<lower=0,upper=1>[N_ind] s_kk;//sensitivity of 1x1 KK
}

transformed data{
   vector[N_ind] inf_eps_a_t;//epg infection intensity of each
  vector<lower=0,upper=1>[N_ind] s_kk_t;//sensitivity of 1x1 KK
  s_kk_t=(s_kk+0.01)/1.02;
   inf_eps_a_t=(inf_eps+0.01)/25;
}

parameters{
 real a0;
 real<lower=0> a1;
 real<lower=0> a2;
 real b0;
 real b1;
}

transformed parameters{
  vector<lower=0,upper=1>[N_ind] mu;
  vector<lower=0>[N_ind] nu; 
  vector<lower=0>[N_ind] alpha;
  vector<lower=0>[N_ind] beta;
  
    mu=inv_logit(a0+a1*exp(log(inf_eps_a_t)./a2));
    for(i in 1:N_ind){
    nu[i]=exp(b0+b1*log(inf_eps_a_t[i]));////worked with log for weird reasons
    }
    alpha=mu ./ nu;
    beta=(1-mu) ./ nu;//replaced nu with 1/nu to not get values in order of 15

}

model{
  a0~normal(-2,4);
  a1~gamma(10,1);
  a2~gamma(5,1);
  b0~normal(0,10);
  b1~normal(0,10);
  s_kk_t~beta(alpha,beta);
}
