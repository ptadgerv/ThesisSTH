data{
 int<lower=1> N_ind; //number of datapoints
 vector<lower=0,upper=1>[N_ind] obs_prev;//epg infection intensity of each dataset
 vector<lower=0,upper=1>[N_ind] true_prev;//sensitivity of 1x1 KK
}


parameters{
 real<lower=0,upper=1> a0;
 real<lower=0> a1;
 real<lower=0> b0;
 real<lower=0> b1;
}

transformed parameters{
  vector<lower=0,upper=1>[N_ind] mu;
  vector<lower=0>[N_ind] nu; 
  vector<lower=0>[N_ind] alpha;
  vector<lower=0>[N_ind] beta;
  
    mu=(2*inv_logit(a1*obs_prev)-1)*(1-a0)+a0;
    for(i in 1:N_ind){
    nu[i]=b0+b1*obs_prev[i];////worked with log for weird reasons
    }
    alpha=mu ./ nu;
    beta=(1-mu) ./ nu;//replaced nu with 1/nu to not get values in order of 15

}

model{
  a0~beta(1,1);
  a1~gamma(1,0.1);
  b0~normal(0,10);
  b1~normal(0,10);
  true_prev~beta(alpha,beta);
}
