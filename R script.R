set.seed(50)
v = 6 ; b = 5 
mu = 0 ; beta = seq(-(b-1)/2,(b-1)/2,length=5); tau = rep(0,v)
sim = 10^5

F_sim = NULL; F_sim_1 = NULL; F_sim_2 = NULL 

for(i in 1:sim)
{
  data_sim = mu + outer(beta,tau,"+") + matrix(rnorm(b*v),nrow = b)
  data_sim[1,1] <- 0
  
  B.1_dash <- sum(data_sim[1,])          # B.1_dash is vector of B1'
  T.1_dash <- sum(data_sim[,1])          # T.1_dash is vector of T1'
  G_dash <- sum(data_sim)                # G_dash is vector of G'
  
  x_hat <- (b*B.1_dash+v*T.1_dash-G_dash)/((b-1)*(v-1))
  x_curl <- B.1_dash/(v-1)
  
  #### approximate test ####
  
  data_sim[1,1] <- x_hat
  B <- rowSums(data_sim)
  T <- colSums(data_sim)
  G <- sum(data_sim)  
  SSBl = sum(B^2)/v - (G^2/(b*v))
  SSTr = sum(T^2)/b - (G^2/(b*v)) ; MSTr = SSTr/(v-1)
  TSS = sum(data_sim^2) - (G^2/(b*v))
  SSE = TSS - SSBl - SSTr ; MSE = SSE/(b*v-b-v)
  F_sim[i] <- MSTr/MSE
  
  sse_hat=SSE
  
  #### More accurate testing procedure-1  ####
  
  SSBl = (B.1_dash^2)/(v-1) + (sum(B[-1]^2)/v) - (G_dash^2/(b*v - 1)) 
  TSS = sum(data_sim^2)-data_sim[1,1]^2 - (G_dash^2/(b*v - 1))
  SSTr = TSS - SSBl - sse_hat ;  MSTr = SSTr/(v-1)
  F_sim_1[i] <- MSTr/(sse_hat/(b*v-b-v))
  
  #### More accurate testing procedure-2  ####
  
  data_sim[1,1] <- x_curl
  B <- rowSums(data_sim)
  T <- colSums(data_sim)
  G <- sum(data_sim)
  
  SSBl = sum(B^2)/v - (G^2/(b*v))
  SSTr = sum(T^2)/b - (G^2/(b*v))
  TSS = sum(data_sim^2) - (G^2/(b*v))
  SSE0 = TSS - SSBl
  F_sim_2[i] <- ((SSE0 - SSE)/(v-1))/(sse_hat/(b*v-b-v))
}

sup = seq(0,6,0.1)
F_ecdf = ecdf(F_sim)
F_ecdf_1 = ecdf(F_sim_1)
F_ecdf_2 = ecdf(F_sim_2)


plot(sup,pf(sup,df1=v-1,df2=b*v-b-v),type = "l",lwd=0.5,lty=2,col="darkgreen",
     xlab="value test statistic, F", ylab="Empirical cdf of the test statistic")
lines(sup,F_ecdf(sup),type = "l",lwd=0.5,lty=2,col="blue")
lines(sup,F_ecdf_2(sup),type = "l",lwd=0.5,col="green")
lines(sup,F_ecdf_1(sup),type = "l",lwd=0.5,lty=2,col="red")
legend(2.5,0.7, legend=c("Approximate test", "More accurate test Procedure-1",
                         "More accurate test Procedure-2"),col=c("blue","green","red"), lty=1:2, cex=0.8)




----------------#####  Power Curves ####--------

set.seed(50)
sim = 10^4; alpha = 0.05 
v0 = 5; b0 = 6
beta0 = seq(-b0,b0,length=b0)
Power_Curve=function(tau=tau0,v=v0,b=b0,mu=0,beta=beta0)
{
  power.vec = NULL; power.vec_1 = NULL; power.vec_2 = NULL
  cut_pt = qf(p = alpha,df1 = (v-1),df2 = (b-1)*(v-1) - 1,lower.tail = FALSE)
  for(i in 1:sim)
  {
    data_sim = mu + outer(beta,tau,"+") + matrix(rnorm(b*v),nrow = b)
    data_sim[1,1] <- 0
    
    B.1_dash <- sum(data_sim[1,])          # B.1_dash is vector of B1'
    T.1_dash <- sum(data_sim[,1])          # T.1_dash is vector of T1'
    G_dash <- sum(data_sim)                # G_dash is vector of G'
    
    x_hat <- (b*B.1_dash+v*T.1_dash-G_dash)/((b-1)*(v-1))
    x_curl <- B.1_dash/(v-1)
    
    #### approximate test ####
    
    data_sim[1,1] <- x_hat
    B <- rowSums(data_sim)
    T <- colSums(data_sim)
    G <- sum(data_sim)  
    SSBl = sum(B^2)/v - (G^2/(b*v))
    SSTr = sum(T^2)/b - (G^2/(b*v)) ; MSTr = SSTr/(v-1)
    TSS = sum(data_sim^2) - (G^2/(b*v))
    SSE = TSS - SSBl - SSTr ; MSE = SSE/(b*v-b-v)
    power.vec[i] = (MSTr/MSE > cut_pt)
    
    sse_hat=SSE
    
    #### More accurate testing procedure-1  ####
    
    SSBl = (B.1_dash^2)/(v-1) + (sum(B[-1]^2)/v) - (G_dash^2/(b*v - 1)) 
    TSS = sum(data_sim^2)-data_sim[1,1]^2 - (G_dash^2/(b*v - 1))
    SSTr = TSS - SSBl - sse_hat ;  MSTr = SSTr/(v-1)
    power.vec_1[i] = (MSTr/(sse_hat/(b*v-b-v)) > cut_pt)
    
    #### More accurate testing procedure-2  ####
    
    data_sim[1,1] <- x_curl
    B <- rowSums(data_sim)
    T <- colSums(data_sim)
    G <- sum(data_sim)
    
    SSBl = sum(B^2)/v - (G^2/(b*v))
    SSTr = sum(T^2)/b - (G^2/(b*v))
    TSS = sum(data_sim^2) - (G^2/(b*v))
    SSE0 = TSS - SSBl
    power.vec_2[i] = (((SSE0 - SSE)/(v-1))/(sse_hat/(b*v-b-v)) > cut_pt)
  }
  power=matrix(c(power.vec,power.vec_1,power.vec_2),nrow=sim)
}


n = 20
sample.tau = NULL
tau_mat = NULL

for(i in 1:n)
{
  sample.tau = seq(-i/10,i/10,length = v0)
  tau_mat = rbind(tau_mat,sample.tau -  mean(sample.tau))
}

sum_tau_sq = apply(tau_mat,MARGIN = 1,FUN = function(x){return(sum(x^2))})

Power = NULL; Power_1 = NULL; Power_2 = NULL
for(i in 1:n)
{
  mat=Power_Curve(tau = tau_mat[i,])
  Power[i] = colMeans(mat)[1]
  Power_1[i] = colMeans(mat)[2]
  Power_2[i] = colMeans(mat)[3]
}

plot(sum_tau_sq,Power,type = "l",xlab=expression(sum(tau^2)))
lines(sum_tau_sq,Power_1,type = "l",lwd=0.5,lty=2,col="blue")
lines(sum_tau_sq,Power_2,type = "l",lwd=0.5,lty=4,col="green")
legend(4,0.7, legend=c("Approximate test","More accurate test Procedure-1",
                       "More accurate test Procedure-2"),col=c("grey","blue","green","red"),
       lty=c(1,2,4,2), cex=0.8)



