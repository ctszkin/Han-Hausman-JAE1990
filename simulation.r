library(mvtnorm)
library(maxLik)
library(compiler)
library(parallel)

dgp <- function(n=500,T=5,b1=c(0.2,0.2),b2=c(0.8,-0.3),d1,d2){
  b <- cbind(b1,b2)
#   
  x <- rmvnorm(n,c(0,0),sigma=matrix(c(5,1.5,1.5,5),ncol=2) )
  e <- rmvnorm(n,sigma=matrix(c(1,0,0,1),ncol=2) )
#   
  d <- x %*%b + e
#   
  t1 <- ceiling(exp((-d[,1])/1.1)+0.05)
  table(t1)
#   t1
#   
  t2 <-  ceiling(exp((-d[,2])/0.9)+0.08)
  table(t2)
#   
  y1_star <- x%*%b1 + e[,1]
  y2_star <- x%*%b2 + e[,2]
#   
#   d1 <- (sort(rweibull(T,1.1)) + 0.05)
#   d2 <- (sort(rweibull(T,0.9)) + 0.08)
#   
#   d1 <- sort(rnorm(T,2,3))
#   d2 <- sort(rnorm(T,2,3))

  t1 <- sapply(y1_star,function(x) min(which(x < d1)))
  t2 <- sapply(y2_star,function(x) min(which(x < d2)))
  
  t1 <- pmin(t1,T)
  t2 <- pmin(t2,T)
  
  data <- list(x=x,mode=(t2<t1)+0, time=pmin(t1,t2))
  data
}

lik_sub_function_1 <- cmpfun(function(x, a, b, c){
  (1 - pnorm(a + (x - b) * c )) * dnorm(x)
})



lik_sub_function_2 <- cmpfun( function(x){
  e1_t <- x[1]
  e1_t1 <- x[2] 
  e2_t <- x[3]
  e2_t1 <- x[4]
  lambda_t <- x[5]
  mode <- x[6]
  
  if(!mode){
    integrate(lik_sub_function_1, lower = e1_t1, upper = e1_t, a=e2_t, b=e1_t, c=lambda_t)$value
  }else{
    integrate(lik_sub_function_1, lower = e2_t1, upper = e2_t, a=e1_t, b=e2_t, c=1/lambda_t)$value
  }
})

loglik <-cmpfun(function(data,b1,b2,d1,d2,cl){
    
  y1_star <- data$x%*%b1
  y2_star <- data$x%*%b2
  
  lambda  <- c(1,head(diff(d2),-1),1) / c(1,head(diff(d1),-1),1)
  
  
  d2_t <- d2[data$time] 
  d1_t <- d1[data$time]
  
  tmp_index <- data$time-1
  tmp_index[tmp_index==0] <- NA
  
  d2_t1 <- d2[tmp_index] 
  d1_t1 <- d1[tmp_index]
  
  d2_t1[is.na(d2_t1)]  <- -100
  d1_t1[is.na(d1_t1)]  <- -100
  
  e2_t <- d2_t - y2_star
  e2_t1 <- d2_t1 -y2_star
  
  e1_t <- d1_t - y1_star
  e1_t1 <- d1_t1 - y1_star
  
  
  lambda_t  <- lambda[data$time]
  
  tmp_data <- cbind(e1_t,e1_t1,e2_t,e2_t1,lambda_t,data$mode)
  
  
  # tmp_data <- subset(tmp_data,subset=data$time!=1)
  
  
  if (missing(cl)){
    sum(log(
      apply(tmp_data,1,lik_sub_function_2)
    ))
  }else{
    sum(log(
      parApply(cl, tmp_data,1,lik_sub_function_2)
    ))
  }
            
            
})

loglik2 <- cmpfun(function(x, data,cl){
  b1 <- x[1:2]
  b2 <- x[3:4]
  d1 <- cumsum(x[5:9])
  d2 <- cumsum(x[10:14])
  loglik(data,b1,b2,d1,d2,cl)
})


cl <- makeCluster(7)
clusterExport(cl,c("lik_sub_function_2","lik_sub_function_1"))



T <- 5
n <- 1500
b1 <- c(0.2,0.2)
b2 <- c(0.8,-0.3)

# tmp <- replicate (10000,sort(rnorm(T,2,3)) )
# d1 <- rowMeans(tmp)
# d2 <- rowMeans(tmp)
#   d1 <- sort(rnorm(T,2,3))
#   d2 <- sort(rnorm(T,2,3))

d1 <- c(-2.7548169, -0.5625119,  0.2822059,  0.4509177,  3.5232749)
d2 <- c(-1.90080534, -0.03773230,  0.08249932,  1.64210103,  3.03850519)
data <- dgp(n=n,T=T,b1=b1,b2=b2,d1,d2)
table(data$time)

start <- c(b1,b2,c(d1[1],diff(d1)),c(d2[1],diff(d2)))

loglik2(start,data,cl)
system.time({
out <- maxLik(loglik2,start=start,method="BFGS", data=data,cl=cl)
})
out
out$estimate[1:4]

# benchmark(loglik2_test1(c(b1,b2,d1,d2),data),loglik2(c(b1,b2,d1,d2),data,cl))



# adding contraint to it
# Ax + B >= 0 
# x is the parameter
# first four is unconstrainted
ineqA <- diag(14)[c(6:9,11:14),]
ineqB <- 0
out2 <- maxLik(loglik2,start=start,method="BFGS", data=data,cl=cl,constraints=list(ineqA=ineqA,ineqB=ineqB))

out2$estimate[1:4]


m <- 500
out<-vector("list",m)
for (i in 53:m){
  data <- dgp(n=n,T=T,b1=b1,b2=b2,d1,d2)
  cat(i,"\n")
  start_time <- proc.time()[3]
  
  tryCatch({
    out[[i]] <- maxLik(loglik2,start=start,method="BFGS", data=data,cl=cl)
  }
  , error = function(e) {}
  )
  
  cat(proc.time()[3]-start_time,"\n")
}
q <- sapply(out[1:282], function(x) x$estimate[1:4])
q <- do.call(rbind,q)

(sqrt(diag(var(q))))
data.frame(mean=colMeans(q),true=c(0.2,0.2,0.8,-0.3),bias=colMeans(q)-c(0.2,0.2,0.8,-0.3),sd=(sqrt(diag(var(q)))))

