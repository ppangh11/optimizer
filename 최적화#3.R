#square loss function (bisection)

square.ftn = function(x,mu){
  fx = mean((x-mu)^2)
  gx = mean(-2*(x-mu))
  return(list(loss=fx,grad=gx))
}

mpg = mtcars$mpg
a = min(mpg)
b = max(mpg)
eps = 1e-7
iter.max = 1000

#square.ftn(mpg,a)
#square.ftn(mpg,b)

square.bi <- c()
for(i in 1:iter.max){
  mid = (a+b)/2
  square.bi = mid
  if(square.ftn(mpg,mid)$grad*square.ftn(mpg,a)$grad>0){
    a = mid
    b = b
  } else{
    a = a
    b = mid
  }
  if(b-a<eps) break
}



#huber loss function(bisection)


huber.ftn = function(x,mu,delta=1){
  fx = mean((x-mu)^2*(abs(x-mu)<delta)+(2*delta*abs(x-mu)-delta^2)*(abs(x-mu)>delta))
  gx = mean(-2*(x-mu)*(abs(x-mu)<=delta)+2*delta*sign(x-mu)*(abs(x-mu)>delta))
  return(list(loss=fx,grad=gx))
}

mpg = mtcars$mpg
a = min(mpg)
b = max(mpg)
eps = 1e-7
iter.max = 1000


huber.bi <- c()
for(i in 1:iter.max){
  mid = (a+b)/2
  huber.bi = mid
  if(huber.ftn(mpg,mid)$grad*huber.ftn(mpg,a)$grad>0){
    a = mid
    b = b
  } else{
    a = a
    b = mid
  }
  if(b-a<eps) break
}

huber.bi


## quantile

quan.ftn = function(alp,mu=0,sigma=1,eps=1e-7,iter.max=1000){
  a = mu-5*sigma
  b = mu+5*sigma
  q.trace = c()
  for(i in 1:iter.max){
    mid = (a+b)/2
    q.trace = c(q.trace,mid)
    if((pnorm(mid,mean=mu,sd=sigma)-alp)*(pnorm(a,mean=mu,sd=sigma)-alp)>0){
      a = mid
      b = b
    }else {
      a = a
      b = mid
    }
    if(b-a<eps) break
  }
  return(q.trace)
}

alp = 0.95
quan.ftn(0.95,0,1)
