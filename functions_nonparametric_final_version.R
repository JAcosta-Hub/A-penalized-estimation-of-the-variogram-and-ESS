## Approximation of maximum distance (observed diameter)
f_hmax <- function(coords){
  rx=range(coords[,1])
  ry=range(coords[,2])
  a=c(rx[1], ry[1])
  b=c(rx[2], ry[2])
  return(as.numeric(dist(rbind(a,b))))
}

# Zeros of Bessel (J0)
tj.opt <- function(m=100){
  tj.aux=c()
  for(j in seq(2,m,by=2) ){
    tj.aux=c(tj.aux, NR(x0=j, fn=besselJ, nu=0)) }
  tj.aux=tj.aux[order(tj.aux)]
  tj.opt=unique(round(tj.aux,6))
  return(tj.opt)
}

# Newton-raphson
NR <- function(x0, fn, gn=NULL, eps=10^(-7), iter.max=1000, is.besselJ0=TRUE,...){
  k=0
  f0=fn(x0,...)
  while(abs(f0)>eps|k<iter.max){
    if(is.besselJ0){g0=-besselJ(x0,nu=1)}else{
      if(is.null(gn)){g0=derivada(x0=x0, fn=fn, eps=eps,...)}else{g0=gn(x0,...)} }
    x1=x0-f0/g0
    if(x1<0){x1=x0-10^(-2)*f0/g0} 
    if(abs(x1-x0)<eps) break
    x0=x1
    f0=fn(x0,...)
    k=k+1
  }
  return(x0)
}

# Determination of nodes in Bessel basis 
fun_tj <- function(m1=10, type="regular"){ 
  tj_opt=tj.opt(m=500) # zeros de bessel
  if(m1>0){
    if(type=="random"){ tj=runif(m1,2,ceiling(tj_opt[m1])); tj=tj[order(tj)] }
    if(type=="regular"){ tj=tj_opt[1:m1] } }else{ tj=NULL } 
  return(tj)
}

# Determination of nodes in exponential basis
fun_uj <- function(m2=10, type="regular", dep.strong=FALSE, plot.it=FALSE){ # uj
  if(m2>0){ 
    if(type=="random"){ uj=runif(m2, 1/60, 1/3); uj=uj[order(uj, decreasing=TRUE)] }
    if(type=="regular"){ uj=seq(1, 1/60, l=m2)/3 } }else{ uj=NULL }
  if(dep.strong){uj=uj*3}
  if(!is.null(uj)&plot.it){
    h=seq(0,5,l=10001)
    plot(NA, xlim=c(0,2), ylim=c(0,1), ylab="")
    for(j in 1:length(uj)){
      points(h, exp(-h/uj[j]), type="l", col=j)
    }
  }
  return(uj)
}

## Calculation of matrix A
Apmatrix <- function(hi, tj, uj){
  if(is.null(tj)){ A1=1 }else{
    A1=cbind(1, -besselJ( outer(hi, tj, FUN='*'), nu=0) ) }
  A2=-exp(-outer(hi, uj, FUN='/'))
  A=cbind(A1,A2)
  return(A)
}

## Calculation of matrix B
Bpmatrix <- function(tj=NULL, uj=NULL, lower=0, upper=1){
  omega2=function(x, a){ a^2*(besselJ(a*x,nu=2)-besselJ(a*x,nu=0))/2  }
  exp2=function(x, a){ exp(-x/a)/a^2}
  
  integrand110 <- function(x, ti){ omega2(x, a=ti)^2 }  
  integrand11 <- function(x, ti, tj){ omega2(x, a=ti)*omega2(x, a=tj) }
  integrand12 <- function(x, ti, uj){ omega2(x, a=ti)*exp2(x, a=uj) }
  integrand22 <- function(x, ui, uj){ exp2(x, a=ui)*exp2(x, a=uj) }
  
  ##-------------------------------##
  ## All block (Bessel-exponential basis)
  if(!is.null(uj)&!is.null(tj)){
    # Block B11
    B11=matrix(0, nc=length(tj)+1, nr=length(tj)+1)
    for(i in 1:(length(tj)-1)){
      for(j in (i+1):length(tj)){
          B11[i+1,j+1]=integrate(integrand11, lower=lower, upper=upper, ti=tj[i], tj=tj[j], subdivisions=10000L)$value }
    }
    B11=t(B11)+B11
    for(i in 1:length(tj)){
      B11[i+1,i+1]=integrate(integrand110, lower=lower, upper=upper, ti=tj[i], subdivisions=10000L)$value }
    
    # Block B12
    B12=matrix(0, nc=length(uj), nr=length(tj)+1)
    for(i in 1:length(tj)){
      for(j in 1:length(uj)){
        B12[i+1,j]=integrate(integrand12, lower=lower, upper=upper, ti=tj[i], uj=uj[j], subdivisions=10000L)$value
      }
    }
    # Block B22  
    aux1=outer(1/uj, 1/uj, FUN='+')
    aux2=outer(uj, uj, FUN='*')
    B22=(1-exp(-aux1))/(aux2^(2)*aux1)
    
    B=rbind(cbind(B11,B12), cbind(t(B12),B22))
  }
  
  ##-------------------------------##
  ## Only Block 11 (Bessel basis)
  if(is.null(uj)&!is.null(tj)){
    B11=matrix(0, nc=length(tj)+1, nr=length(tj)+1)
    for(i in 1:(length(tj)-1)){
      for(j in (i+1):length(tj)){
          B11[i+1,j+1]=integrate(integrand11, lower=lower, upper=upper, ti=tj[i], tj=tj[j], subdivisions=1000000L)$value }
    }
    B11=t(B11)+B11
    for(i in 1:length(tj)){
      B11[i+1,i+1]=integrate(integrand110, lower=lower, upper=upper, ti=tj[i], subdivisions=1000000L)$value }
    B=B11
  }
  
  ##-------------------------------##
  ## Only Block22 (exponential basis)
  if(!is.null(uj)&is.null(tj)){
    aux1=outer(1/uj, 1/uj, FUN='+')
    aux2=outer(uj, uj, FUN='*')
    B22=(1-exp(-aux1))/(aux2^(2)*aux1)
    
    B=rbind(0, cbind(0,B22))
  }
  if(is.null(uj)&is.null(tj)){ # Both null (pure nugget)
    B=0 }
  return(B)
}


## Objetive function for estimated Lambda
lambda.opt <- function(lambda, m1=4, m2=4, tj=NULL, uj=NULL, Ap=NULL, Bp=NULL, hi=NULL, y=NULL, vario, hmax){
  require(pracma)
  if(is.null(hi)){hi=vario$u}
  if(is.null(y)){y=vario$v}
  if(is.null(Ap)){
    if(is.null(tj)){ tj=fun_tj(m1=m1) }else{ m1=length(tj) } 
    if(is.null(uj)){ uj=fun_uj(m2=m2) }else{ m2=length(uj) } 
    Ap=Apmatrix(tj=tj, uj=uj, hi=hi/hmax) }
  if(is.null(Bp)){
    if(is.null(tj)){ tj=fun_tj(m1=m1) }else{ m1=length(tj) } 
    if(is.null(uj)){ uj=fun_uj(m2=m2) }else{ m2=length(uj) } 
    Bp=Bpmatrix(tj=tj, uj=uj)/hmax^3
  }  
  nm=dim(Ap)
  S0=Ap%*%pracma::pinv(t(Ap)%*%Ap+lambda*Bp)%*%t(Ap)
  df0=sum(diag(S0)); #df
  GCV=mean(((diag(nm[1])-S0)%*%y)^2/(1-mean(diag(S0)))^2)
  return(GCV)
}

## Objetive function for estimated variogram 
Qp_v2 <- function(z, y, A, B, lambda, wi){
  aux1=A%*%z
  aux2=B%*%z
  res = sum(wi*(y-aux1)^2)+lambda*sum(z*aux2) 
  return(res)
}

## Objetive function for estimated variogram: Force nugget to 0 or fixed
Qp0_v2 <- function(z, y, A, B, lambda, wi){
  A0=A[,-1]
  B0=B[-1,-1]
  aux1=(A0+1)%*%z
  aux2=B0%*%z
  res = sum(wi*(y-aux1)^2)+lambda*sum(z*aux2)
  return(res)
}

## Estimation Function
NonParametric2 <- function(datos, coords=NULL, lambda=0, tj=NULL, uj=NULL, 
                           m=10, m2=10, hmax=NULL, A=NULL, B=NULL, type="regular", type.weigth=1,
                           fixnugget=FALSE, nugget=0){
  
  if(is.null(hmax)&is.null(coords)){print("ingrese hmax o las coordenadas")}
  if(is.null(hmax)&!is.null(coords)){hmax=f_hmax(coords)}

  ## format of data
  if(class(datos)=="variogram"){ 
    vg=datos 
    y=vg$v
    hi=vg$u
    vg.n=vg$n
    vg.sd=vg$sd }
  if(class(datos)=="geodata"){ 
    vg=geoR::variog(datos, max.dist=hmax/3 )
    y=vg$v
    hi=vg$u
    vg.n=vg$n
    vg.sd=vg$sd
  }
  if(class(datos)=="GeoVariogram"){
    y=datos$variograms 
    hi=datos$centers
    vg.n=datos$lenbins}
  if(class(datos)=="data.frame"&!is.null(coords)){
    vg=GeoModels::GeoVariogram(data=datos,coordx=coords, maxdist=hmax/3, numbins=13)
    y=vg$variograms 
    hi=vg$centers
    vg.n=vg$lenbins
  }
  
  # weights
  if(type.weigth==1){ wi=1 }
  if(type.weigth==2){ wi=(1/vg.n)/sum(1/vg.n) }
  if(type.weigth==3){ wi=(1/vg.sd)/sum(1/vg.sd) }
  if(type.weigth==4){ wi=(1/vg.sd^2)/sum(1/vg.sd^2) }
  if(type.weigth==5){ wi=(1/((vg.sd/y)*vg.n))/sum(1/((vg.sd/y)*vg.n)) }
  
  # nodes
  if(is.null(tj)){ tj=fun_tj(m1=m, type=type) }else{ m=length(tj) } 
  if(is.null(uj)){ uj=fun_uj(m2=m2, type=type) }else{ m2=length(uj) } 
  
  # Calculations of matrix A and B
  if(is.null(A)){A=Apmatrix(tj=tj, uj=uj, hi=hi/hmax)}
  if(is.null(B)){B=Bpmatrix(tj=tj, uj=uj)/hmax^3}
  
  # initial value
  z0=c(max(y), rep(0,m), rep(0,m2))

  # Estimations
  if(!fixnugget){ # estimation with non-fixed nugget
      opt=optim(par=z0, fn=Qp_v2, lower=rep(0,m+m2+1), method="L-BFGS-B", 
                y=y, A=A, B=B, lambda=lambda, wi=wi)
      z1=opt$par}else{ z1=z0+c(0, rep(abs(max(y)-nugget),m+m2)/(m+m2) ) }
    
  if(z1[1]<sum(z1[-1])|fixnugget){ # estimation with fixed nugget
      print("nugget is forced to 0 or fixed value")
      z0=z1[-1]
      y=y-nugget
      opt=optim(par=z0, fn=Qp0_v2, lower=rep(0,m+m2), method="L-BFGS-B", 
                y=y, A=A, B=B, lambda=lambda, wi=wi)
      z1=opt$par
      opt$par <- c(sum(z1)+nugget, z1)
    }
  
  return(list(opt=opt, tj=tj, uj=uj))
}



## simulation result (in parallel) with optimum lambda
res3.sim <- function(vario.datos.sim, hmax, n.cores=NULL, 
                     lambda=NULL, m.sim=NULL, m1=20, m2=20, 
                     lambda.ini=0.01, fixnugget=FALSE){
  start_time <- Sys.time()
  # initialization
  if(is.null(m.sim)){m.sim=length(vario.datos.sim)}
  if(is.null(n.cores)){n.cores=detectCores()}
  print(paste("the number of simulations is", m.sim))
  print(paste("number of cores used", n.cores))
  jvalues=1:m.sim
  
  Bp1=Bpmatrix(tj=fun_tj(m1=m1), uj=fun_uj(m2=0))/hmax^3
  Bp2=Bpmatrix(tj=fun_tj(m1=0), uj=fun_uj(m2=m2))/hmax^3
  Bp3=Bpmatrix(tj=fun_tj(m1=ceiling(m1/2)), uj=fun_uj(m2=ceiling(m2/2)))/hmax^3
  print("the dimension of B is"); print(dim(Bp1))
  print(paste("The nugget is estimated:", !fixnugget))
  
  # Creation of the cluster
  cl <- makeCluster(n.cores)
  clusterExport(cl, c('hmax', 'jvalues', 'm1', 'm2', 'lambda.opt', 
                      'tj.opt', 'NR', 'fun_tj', 'fun_uj', #'lambda.seq', 'l',
                      'Apmatrix', 'NonParametric2', 'Qp_v2', 'Qp0_v2'))
  res.lambda <- parLapply(cl = cl, # Cluster
                          jvalues, 
                          function(j,...){
                            hi=vario.datos.sim[[j]]$u
                            yi=vario.datos.sim[[j]]$v
                            
                            Ap1=Apmatrix(hi=hi/hmax, tj=fun_tj(m1=m1), uj=fun_uj(m2=0))
                            Ap2=Apmatrix(hi=hi/hmax, tj=fun_tj(m1=0), uj=fun_uj(m2=m2))
                            Ap3=Apmatrix(hi=hi/hmax, tj=fun_tj(m1=ceiling(m1/2)), uj=fun_uj(m2=ceiling(m2/2)))
                            
                            if(is.null(lambda)){
                              aux1_select <- optim(par=lambda.ini, fn=lambda.opt, method="L-BFGS-B", 
                                                   y=yi, hi=hi, hmax=hmax, lower=10^(-3), upper=Inf, Ap=Ap1, Bp=Bp1)
                              
                              aux2_select <- optim(par=lambda.ini, fn=lambda.opt, method="L-BFGS-B", 
                                                   y=yi, hi=hi, hmax=hmax, lower=10^(-3), upper=Inf, Ap=Ap2, Bp=Bp2)
                              
                              aux3_select <- optim(par=lambda.ini, fn=lambda.opt, method="L-BFGS-B", 
                                                   y=yi, hi=hi, hmax=hmax, lower=10^(-3), upper=Inf, Ap=Ap3, Bp=Bp3)
                              
                              
                              res1.aux=NonParametric2(datos=vario.datos.sim[[j]], lambda=aux1_select$par, 
                                                      m=m1, m2=0, hmax=hmax, A=Ap1, B=Bp1, fixnugget=fixnugget)
                              
                              res2.aux=NonParametric2(datos=vario.datos.sim[[j]], lambda=aux2_select$par, m=0, 
                                                      m2=m2, hmax=hmax, A=Ap2, B=Bp2, fixnugget=fixnugget)
                              
                              res3.aux=NonParametric2(datos=vario.datos.sim[[j]],  lambda=aux3_select$par, 
                                                      m=ceiling(m1/2), m2=ceiling(m2/2), A=Ap3, B=Bp3,
                                                      fixnugget=fixnugget)
                              
                              return(list(res.estim=list(bessel=res1.aux, expon=res2.aux, ambos=res3.aux),
                                          opt.lambda=rbind(bessel=c(aux1_select$par, aux1_select$value), 
                                                           expon=c(aux2_select$par, aux2_select$value), 
                                                           ambos=c(aux3_select$par, aux3_select$value) )))
                              
                            }else{
                              res1.aux=NonParametric2(datos=vario.datos.sim[[j]], lambda=lambda, m=m1, m2=0, hmax=hmax, A=Ap1, B=Bp1, fixnugget=fixnugget)
                              res2.aux=NonParametric2(datos=vario.datos.sim[[j]], lambda=lambda, m=0, m2=m2, hmax=hmax, A=Ap2, B=Bp2, fixnugget=fixnugget)
                              res3.aux=NonParametric2(datos=vario.datos.sim[[j]], lambda=lambda, hmax=hmax, A=Ap3, B=Bp3, 
                                                      m=ceiling(m1/2), m2=ceiling(m2/2), fixnugget=fixnugget)

                              return(list(res.estim=list(bessel=res1.aux, expon=res2.aux, ambos=res3.aux),
                                          opt.lambda=rbind(bessel=lambda, 
                                                           expon=lambda, 
                                                           ambos=lambda )))
                            }
                            
                          })
  
  stopCluster(cl) # close the clusters
  finish_time <- Sys.time()
  tmp = finish_time - start_time
  print(tmp) # Time of processing
  return(list(res=res.lambda, time=tmp))
}

## Auxiliary function for ML estimation
ifmax<-function(a){
  aux.loglik=c()
  if(length(a)==1){ return(a[[1]])}else{
    for(i in 1:length(a)){ 
      if(inherits(a[[i]], "try-error")){ aux.loglik=c(aux.loglik, NA); next }
      aux.loglik=c(aux.loglik, a[[i]]$loglik) }
    res=a[[which.max(aux.loglik)]]
    return(res)
  }
}


## simulation result (in parallel) for Maximum Likelihood
res2_ML.sim <- function(geo.datos.sim, n.cores=NULL, m.sim=NULL, fix.nugget=TRUE,
                        cov.model="exponential", nugget=0.01, ini.cov.pars, kappa=0.5,...){
  # initialization
  if(is.null(m.sim)){m.sim=length(geo.datos.sim)}
  if(is.null(n.cores)){n.cores=detectCores()}
  print(paste("the number of simulations is", m.sim))
  print(paste("number of cores used", n.cores))
  jvalues=1:m.sim
  
  # Creation of the cluster
  cl <- makeCluster(n.cores)
  clusterExport(cl, c('jvalues', 'ifmax'))
  start_time <- Sys.time()
  
  res.ml<- parLapply(cl = cl, # Cluster
                     jvalues,
                     function(j,...){
                       res.aux.ML=list()
                       res.aux.REML=list()
                       for(i in 1:length(cov.model)){
                         res.aux.ML[[i]]=try(geoR::likfit(geo.datos.sim[[j]], nugget=nugget, cov.model=cov.model[i], fix.nugget=fix.nugget,
                                                          ini.cov.pars=ini.cov.pars, lik.method="ML", kappa=kappa, messages=FALSE,...),
                                             silent = TRUE)
                         res.aux.REML[[i]]=try(geoR::likfit(geo.datos.sim[[j]], nugget=nugget, cov.model=cov.model[i], fix.nugget=fix.nugget,
                                                            ini.cov.pars=ini.cov.pars, lik.method="REML", kappa=kappa, messages=FALSE,...),
                                               silent = TRUE)
                       }
                       return( list( ML = ifmax(res.aux.ML), REML=ifmax(res.aux.REML) ) )
                     })
  finish_time <- Sys.time()
  stopCluster(cl) # close the clusters
  
  tmp = finish_time - start_time
  print(tmp) # Time of processing
  return(list(res=res.ml, time=tmp))
}



## Auxiliary function for OLS estimation
ifmax2<-function(a){
  aux.ss=c()
  if(length(a)==1){ return(a[[1]])}else{
    for(i in 1:length(a)){ 
      if(inherits(a[[i]], "try-error")){ aux.ss=c(aux.ss, NA); next }
      aux.ss=c(aux.ss, summary(a[[i]])$sum.of.squares) }
    res=a[[which.max(aux.ss)]]
    return(res)
  }
}

## simulation result (in parallel) for Least Square
res2_OLS.sim <- function(vario.datos.sim, n.cores=NULL, m.sim=NULL, 
                         cov.model="exponential", nugget=0.01, ini.cov.pars, kappa=0.5,...){
  # initialization
  if(is.null(m.sim)){m.sim=length(vario.datos.sim)}
  if(is.null(n.cores)){n.cores=detectCores()}
  print(paste("the number of simulations is", m.sim))
  print(paste("number of cores used", n.cores))
  jvalues=1:m.sim
  
  # Creation of the cluster
  cl <- makeCluster(n.cores)
  clusterExport(cl, c('jvalues', 'ifmax2'))
  start_time <- Sys.time()
  
  res.variofit <- parLapply(cl = cl, # Cluster
                            jvalues, 
                            function(j,...){
                              res.aux.OLS=list()
                              res.aux.WLS=list()
                              for(i in 1:length(cov.model)){
                                res.aux.OLS[[i]]=try(geoR::variofit(vario.datos.sim[[j]], nugget=nugget, cov.model=cov.model[i], 
                                                                    ini.cov.pars=ini.cov.pars, weights="equal", kappa=kappa, messages=FALSE,...),
                                                     silent = TRUE)
                                res.aux.WLS[[i]]=try(geoR::variofit(vario.datos.sim[[j]], nugget=nugget, cov.model=cov.model[i], 
                                                                    ini.cov.pars=ini.cov.pars, weights="npairs", kappa=kappa, messages=FALSE,...),
                                                     silent = TRUE)
                              }
                              return( list( OLS = ifmax2(res.aux.OLS), WLS=ifmax2(res.aux.WLS) ) )
                            })
  finish_time <- Sys.time()
  stopCluster(cl) # close the clusters
  
  tmp = finish_time - start_time
  print(tmp) # Time of processing
  return(list(res=res.variofit, time=tmp))
}


## Result curve
line_result2 <- function(z, tj, uj, hmax=1, lag.max=NULL, print.cond=FALSE, h=NULL){
  if( is.null(lag.max) ){ lag.max=hmax }
  if(print.cond){
    print(paste("Fulfills condition:",z[1]>=sum(z[-1]), "nugget =", z[1]-sum(z[-1]) ) ) 
  }
  y2=c()
  if( is.null(h) ){ h=seq(0,lag.max,l=1001) }
  for(i in 1:length(h)){
    y2[i] = z[1]-sum(z[2:(length(tj)+1)]*besselJ(tj*h[i]/hmax,nu=0))-
      sum(z[-(1:(length(tj)+1))]*exp(-h[i]/(uj*hmax)))
  }
  return(data.frame(h=h,gamma=y2))
}

##------------------------------------
## Realizations graph
plot1.ejem <- function(j=1, vario1.datos.sim, vario2.datos.sim, vario3.datos.sim,
                       res.lambda1, res.lambda2, res.cloud, res.lik, res.vario1, res.vario2, fcov.sim, hmax=0.5){
  h=seq(0.0001, hmax, by=0.01)
  tipo_funcion=paste(c("Bessel", "Exponential", "Bessel-Expon"),"basis", sep=" ")
  
  for(k in 1:3){
    plot(vario2.datos.sim[[j]], pch=20, col="gray60", ylim=c(0, 1.75), xlim=c(0,hmax), cex=0.25, main=tipo_funcion[k],
         ylab=bquote(gamma(h)), xlab=bquote(h))
    points(vario1.datos.sim[[j]]$u, vario1.datos.sim[[j]]$v, pch=20, col=1, cex=0.75)
    points(vario3.datos.sim[[j]]$u, vario3.datos.sim[[j]]$v, pch=20, col=2, cex=0.75)
    
    # NP-cloud
    lines(line_result2(z=res.cloud$res[[j]]$res.estim[[k]]$opt$par, 
                       tj=res.cloud$res[[j]]$res.estim[[k]]$tj,
                       uj=res.cloud$res[[j]]$res.estim[[k]]$uj, hmax=hmax, h=h), col="gray40", lwd=2.0, lty=2)
    
    lines(h, fcov.sim(0)-fcov.sim(h), col=6, lwd=1.5) # TRUE
    
    # MV-OLS
    lines(res.lik$res[[j]]$ML, col=3, lwd=2.0, max.dist=hmax, lty=2)
    lines(res.vario1$res[[j]]$OLS, col=4, lwd=2.0, max.dist=hmax, lty=2)
    lines(res.vario2$res[[j]]$OLS, col=5, lwd=2.0, max.dist=hmax, lty=2)
    
    # NP
    lines(line_result2(z=res.lambda1$res[[j]]$res.estim[[k]]$opt$par, 
                       tj=res.lambda1$res[[j]]$res.estim[[k]]$tj,
                       uj=res.lambda1$res[[j]]$res.estim[[k]]$uj, hmax=hmax, h=h), col=1, lwd=2.0, lty=2) 
    lines(line_result2(z=res.lambda2$res[[j]]$res.estim[[k]]$opt$par, 
                       tj=res.lambda2$res[[j]]$res.estim[[k]]$tj,
                       uj=res.lambda2$res[[j]]$res.estim[[k]]$uj, hmax=hmax, h=h), col=2, lwd=2.0, lty=2) 
    
    legend("bottomright", legend=c("True","ML", "OLS-bins48", "OLS-bins16", "NP-bins48", "NP-bins16", "NP-Cloud"), 
           lty=c(1,2,2,2,2,2,2), col=c(6,3:5,1,2,"gray40"), lwd=c(1.5,2,2,2,2,2,2,2))
  }
} 

plot2.ejem <- function(j=1, vario1.datos.sim, vario2.datos.sim, vario3.datos.sim,
                       res.lambda1, res.lambda2, res.cloud, fcov.sim, hmax=1,...){
  h=seq(0.0001, hmax, by=0.01)
  tipo_funcion=paste(c("Bessel", "Exponential", "Bessel-Expon"),"basis", sep=" ")
  
  for(k in 1:3){
    plot(vario2.datos.sim[[j]], pch=20, col="gray60", ylim=c(0, 1.85), xlim=c(0,hmax), cex=0.25, main=tipo_funcion[k],
         ylab=bquote(gamma(h)), xlab=bquote(h))
    points(vario1.datos.sim[[j]]$u, vario1.datos.sim[[j]]$v, pch=20, col=1, cex=0.75)
    points(vario3.datos.sim[[j]]$u, vario3.datos.sim[[j]]$v, pch=20, col=2, cex=0.75)
    
    #NP-Cloud
    lines(line_result2(z=res.cloud$res[[j]]$res.estim[[k]]$opt$par, 
                       tj=res.cloud$res[[j]]$res.estim[[k]]$tj,
                       uj=res.cloud$res[[j]]$res.estim[[k]]$uj, hmax=hmax, h=h,...), col="gray40", lwd=2.0, lty=1)
    #TRUE
    lines(h, fcov.sim(0)-fcov.sim(h), col=4, lwd=1.5)
    
    # NP
    lines(line_result2(z=res.lambda1$res[[j]]$res.estim[[k]]$opt$par, 
                       tj=res.lambda1$res[[j]]$res.estim[[k]]$tj,
                       uj=res.lambda1$res[[j]]$res.estim[[k]]$uj, hmax=hmax, h=h,...), col=1, lwd=2.0, lty=1) 
    lines(line_result2(z=res.lambda2$res[[j]]$res.estim[[k]]$opt$par, 
                       tj=res.lambda2$res[[j]]$res.estim[[k]]$tj,
                       uj=res.lambda2$res[[j]]$res.estim[[k]]$uj, hmax=hmax, h=h,...), col=2, lwd=2.0, lty=1) 
    
    legend("bottomright", legend=c("True", "NP-bins48", "NP-bins16", "NP-Cloud"), 
           lty=1, col=c(4,1,2,"gray40"), lwd=c(1.5,2,2,2))
  }
} 



## Performance metrics
metrics3.sim <- function(res.sim, hmax=1, fcov.sim, type="lambda", h=NULL, ...){
  n.sim=length(res.sim)
  y_real=sapply(h, fvario,  fcov=fcov.sim)

  if(type=="lambda"){
    n.tipo=length(res.sim[[1]]$res.estim)
#    print(cbind(n.sim=n.sim, n.tipo=n.tipo, n.h=length(h), hmax=hmax))
    
    res_rmse <- res_mae <- array(NA, dim=c(n.sim, n.tipo) )
    res_rmse_estand <- res_mae_estand <- array(NA, dim=c(n.sim, n.tipo) )
    for(k in 1:n.tipo){
      for(j in 1:n.sim){
        y_hat=line_result2(z=res.sim[[j]]$res.estim[[k]]$opt$par, 
                           tj=res.sim[[j]]$res.estim[[k]]$tj,
                           uj=res.sim[[j]]$res.estim[[k]]$uj, hmax=hmax, h=h,...)$gamma
        
        res_rmse[j,k] <- sqrt(mean((y_real-y_hat)^2))
        res_mae[j,k] <- mean(abs(y_real-y_hat))
        res_rmse_estand[j,k] <- sqrt(mean(((y_real-y_hat)/y_real)^2))
        res_mae_estand[j,k] <- mean(abs((y_real-y_hat)/y_real))
      }
      
    }
    
  }
  if(type=="cloud"){
    n.tipo=length(res.sim[[1]]$res.estim)
 #   print(cbind(n.sim=n.sim, n.tipo=n.tipo, n.h=length(h), hmax=hmax))
    
    res_rmse <- res_mae <- array(NA, dim=c(n.sim, n.tipo) )
    res_rmse_estand <- res_mae_estand <- array(NA, dim=c(n.sim, n.tipo) )
    for(k in 1:n.tipo){
      for(j in 1:n.sim){
        y_hat=line_result2(z=res.sim[[j]]$res.estim[[k]]$opt$par, 
                           tj=res.sim[[j]]$res.estim[[k]]$tj,
                           uj=res.sim[[j]]$res.estim[[k]]$uj, hmax=hmax, h=h,...)$gamma
        res_rmse[j,k] <- sqrt(mean((y_real-y_hat)^2))
        res_mae[j,k] <- mean(abs(y_real-y_hat))
        res_rmse_estand[j,k] <- sqrt(mean(((y_real-y_hat)/y_real)^2))
        res_mae_estand[j,k] <- mean(abs((y_real-y_hat)/y_real))
      }
      
    }
    
  }
  if(type=="lik"|type=="vario"){
    n.tipo=length(res.sim[[1]])
 #   print(cbind(n.sim=n.sim, n.tipo=n.tipo, n.h=length(h)))
    
    res_rmse <- res_mae <- array(NA, dim=c(n.sim, n.tipo) )
    res_rmse_estand <- res_mae_estand <- array(NA, dim=c(n.sim, n.tipo) )
    for(k in 1:n.tipo){
      for(j in 1:n.sim){
        if(inherits(res.sim[[j]][[k]], "try-error") ){y_hat=NA}else{
          y_hat=res.sim[[j]][[k]]$nugget+res.sim[[j]][[k]]$cov.pars[1]-
            geoR::cov.spatial(h, cov.model = res.sim[[j]][[k]]$cov.model,
                              cov.pars = res.sim[[j]][[k]]$cov.pars,
                              kappa = res.sim[[j]][[k]]$kappa)
        }
        res_rmse[j,k] <- sqrt(mean((y_real-y_hat)^2))
        res_mae[j,k] <- mean(abs(y_real-y_hat))
        res_rmse_estand[j,k] <- sqrt(mean(((y_real-y_hat)/y_real)^2))
        res_mae_estand[j,k] <- mean(abs((y_real-y_hat)/y_real))
      }
    }
    
  }
 
  rmse=apply(res_rmse,2,mean, na.rm=TRUE)
  mae=apply(res_mae,2,mean, na.rm=TRUE)
  rmse_estand=apply(res_rmse_estand,2,mean, na.rm=TRUE)
  mae_estand=apply(res_mae_estand,2,mean, na.rm=TRUE)
  
  if(type=="lambda"|type=="cloud"){
    names(rmse) <- names(mae) <-  c("Bessel", "Expon", "Both")
    names(rmse_estand) <- names(mae_estand) <-  c("Bessel", "Expon", "Both")
  }
  if(type=="lik"){
    names(rmse) <- names(mae) <-  c("ML", "REML")
    names(rmse_estand) <- names(mae_estand) <-  c("ML", "REML")
  }
  if(type=="vario"){
    names(rmse) <- names(mae) <-  c("OLS", "WLS")
    names(rmse_estand) <- names(mae_estand) <-  c("OLS", "WLS")
  }
  
  All=data.frame(RMSE=rmse, MAE=mae, RMSE.Estand=rmse_estand, MAE.Estand=mae_estand)
  return( list( RMSE_ind=res_rmse, MAE_ind=res_mae,
                RMSE_ind_estand=res_rmse_estand, MAE_ind_estand=res_mae_estand,
                RMSE=rmse, MAE=mae, RMSE.Estand=rmse_estand, MAE.Estand=mae_estand, All=All ) )
}



## Integral Errors
IntError <- function(res, hgrid, fcov.sim=fcov0.sim, hmax=hmax,...){
  m=length(res)
  IntError=list()
  for(j in 1:m){
    print(paste('Performance for the model', j, 'of', m))
    IntError[[j]]=metrics3.sim(res.sim=res[[j]]$res, 
                               hmax=hmax, fcov.sim=fcov.sim[[j]], h=hgrid,...)
  }
  return(IntError)
}

## Summary Tables
res_table <- function(IntErr){
  tab_ErrInt1 <- rbind()
  tab_ErrInt2 <- rbind()
  for(j in 1:8){
    aux = cbind(IntErr[['all']][[j]]$All[c("RMSE","MAE")],
                IntErr[['min']][[j]]$All[c("RMSE","MAE")],
                IntErr[['inter']][[j]]$All[c("RMSE","MAE")],
                IntErr[['max']][[j]]$All[c("RMSE","MAE")])[,c(1,3,5,7,2,4,6,8)] 
    if(j==1|j==3|j==5|j==7){ tab_ErrInt1 <- rbind(tab_ErrInt1, aux) }
    if(j==2|j==4|j==6|j==8){ tab_ErrInt2 <- rbind(tab_ErrInt2, aux) }
  }
  
  tab_ErrInt1_estand <- rbind()
  tab_ErrInt2_estand <- rbind()
  for(j in 1:8){
    aux = cbind(IntErr[['all']][[j]]$All[c("RMSE.Estand","MAE.Estand")],
                IntErr[['min']][[j]]$All[c("RMSE.Estand","MAE.Estand")],
                IntErr[['inter']][[j]]$All[c("RMSE.Estand","MAE.Estand")],
                IntErr[['max']][[j]]$All[c("RMSE.Estand","MAE.Estand")])[,c(1,3,5,7,2,4,6,8)] 
    if(j==1|j==3|j==5|j==7){ tab_ErrInt1_estand <- rbind(tab_ErrInt1_estand, aux) }
    if(j==2|j==4|j==6|j==8){ tab_ErrInt2_estand <- rbind(tab_ErrInt2_estand, aux) }
  }
  
  
  return(list(tab_ErrInt1, tab_ErrInt2, tab_ErrInt1_estand, tab_ErrInt2_estand))
}


###-------------------------------------------------------------------
## Functions required for ESS calculation

# function that finds the SILL
C0.fun <- function(res){
  m=length(res$res)
  c0=c()
  for(k in 1:3){
    c01=c()
    for(j in 1:m){
      c01[j]=res$res[[j]]$res.estim[[k]]$opt$par[1]
    }
    c0=cbind(c0, c01)
  }
  colnames(c0) <- c("Bessel", "Exponential", "Both")
  return(c0)
}

# function that finds the nugget
nugget.fun <- function(res){
  m=length(res$res)
  c0=c()
  for(k in 1:3){
    c01=c()
    for(j in 1:m){
      c01[j]=res$res[[j]]$res.estim[[k]]$opt$par[1]-sum(res$res[[j]]$res.estim[[k]]$opt$par[-1])
    }
    c0=cbind(c0, c01)
  }
  colnames(c0) <- c("Bessel", "Exponential", "Both")
  return(c0)
}

# estimated variogram
fit.gamma <- function(h, estim, hmax){
  if(!is.matrix(h)){
    fit.gamma <- line_result2(z=estim$opt$par,
                              tj=estim$tj,
                              uj=estim$uj, h=h, hmax=hmax)$gamma 
  }else{
    dm=dim(h)
    h2=c(h)
    fit.gamma <- line_result2(z=estim$opt$par,
                              tj=estim$tj,
                              uj=estim$uj, h=h2, hmax=hmax)$gamma 
    fit.gamma <- matrix(fit.gamma, nr=dm[1], nc=dm[2])
  }
  return(fit.gamma)
}

## ESS from variogram
ESS_gamma <- function(h, estim, n=NULL, hmax=NULL){
  if(is.null(n)&is.matrix(h)){n=ncol(h)}
  if(is.null(n)&!is.matrix(h)){ print("must enter the sample size n"); stop }
  if(is.null(hmax)){print("must enter the maximum distance from the region"); stop }
  
  c0=estim$opt$par[1]
  nugget=estim$opt$par[1]-sum(estim$opt$par[-1])
  
  gamma_fit=fit.gamma(h=h, estim=estim, hmax=hmax)
  
  if(is.matrix(h)){
    ess=c0/(c0+nugget/n-mean(gamma_fit))
  }else{
    if(class(h)=="dist"){h=as.numeric(h)}
    ess=c0/(c0-2*sum(gamma_fit)/n^2)
  }
  
  return(ess)
}

## Calculation of ESS from simulation estimates
ESS_NP_sim <- function(res.lambda){
  start_time <- Sys.time()
  cl <- makeCluster(4)
  clusterExport(cl, c('ESS_gamma', 'dist2', 'fit.gamma', 'line_result2', 'res.lambda', 'hmax') )
  
  ESS_estim <- parLapply(cl = cl, # Cluster
                         1:8, # vector a recorrer
                         function(i){
                           require(parallel)
                           ESS_1=c()
                           for(k in 1:3){
                             cl2 <- makeCluster(10)
                             clusterExport(cl2, c('ESS_gamma', 'dist2', 'res.lambda', 'fit.gamma', 'line_result2', 'hmax'))
                             
                             ESS_2 <- parSapply(cl = cl2, # Cluster
                                                1:1000, 
                                                function(j){
                                                  return(ESS_gamma(h=dist2, estim=res.lambda[[i]]$res[[j]]$res.estim[[k]], n=200, hmax=hmax) ) })
                             stopCluster(cl2) # cerrar los cluster
                             ESS_1=cbind(ESS_1, ESS_2)
                           }
                           colnames(ESS_1) <- c("Bessel", "Exponential", "Both")
                           return(ESS_1)
                         })
  stopCluster(cl) # cerrar los cluster
  finish_time <- Sys.time()
  (t2 = finish_time - start_time)
  print(t2)
  return(ESS_estim)
}


## verificar desde acá

## Estimaciones de la matriz de correlación en paralelo
Rh2_np <- function(H, model, type="lambda.opt", type.function="bessel", 
                   m.sim=NULL, n.cores=NULL, hmax){
  # parámetros iniciales
  if(is.null(m.sim)){m.sim=length(model)}
  if(is.null(n.cores)){n.cores=detectCores()}
  print(paste("la cantidad de simulaciones es", m.sim))
  print(paste("numero de cores utilizados", n.cores))
  print(paste("La cantidad de datos georreferenciados son", ncol(H)))
  jvalues=1:m.sim
  
  # Creacion del cluster
  cl <- makeCluster(n.cores)
  clusterExport(cl, c('jvalues', 'H', 'line_result2', 'hmax'))
  start_time <- Sys.time()
  
  Rh<- parLapply(cl = cl, # Cluster
                 jvalues, # vector a recorrer
                 function(j,...){
                   if(type=="lambda.opt"){
                     C0=model[[j]]$res.estim[[type.function]]$opt$par[1]
                     nugget <- C0-sum(model[[j]]$res.estim[[type.function]]$opt$par[-1])
                     gamma1=array(NA, dim=dim(H))
                     for(i in 1:ncol(H)){
                       gamma1[i,] <- line_result2(z=model[[j]]$res.estim[[type.function]]$opt$par, 
                                                  tj=model[[j]]$res.estim[[type.function]]$tj,
                                                  uj=model[[j]]$res.estim[[type.function]]$uj, h=H[i,], hmax=hmax)$gamma
                       gamma1[i, gamma1[i,]>C0+nugget] <- C0+nugget
                       #                       gamma1[i, H[i,] > hmax]  <- C0+nugget
                     }
                     Ch <- C0-gamma1+nugget*diag(nrow(H))
                     Rh <- Ch/Ch[1,1]
                     Rh[Rh<0]<-0
                   }
                   if(type=="cloud"){
                     C0=model[[j]][[type.function]]$opt$par[1]
                     nugget <- line_result2(z=model[[j]][[type.function]]$opt$par, 
                                            tj=model[[j]][[type.function]]$tj,
                                            uj=model[[j]][[type.function]]$uj, h=0, hmax=hmax)$gamma
                     gamma1=array(NA, dim=dim(H))
                     for(i in 1:ncol(H)){
                       gamma1[i,] <- line_result2(z=model[[j]][[type.function]]$opt$par, 
                                                  tj=model[[j]][[type.function]]$tj,
                                                  uj=model[[j]][[type.function]]$uj, h=H[i,], hmax=hmax)$gamma
                       gamma1[i, gamma1[i,] > C0+nugget] <- C0+nugget
                       #                      gamma1[i, H[i,] > hmax]  <- C0+nugget
                     }
                     Ch <- C0-gamma1+nugget*diag(nrow(H))
                     Rh <- Ch/Ch[1,1]
                     Rh[Rh<0]<-0
                   }
                   if(type=="loglik"){
                     if(inherits(model[[j]][[type.function]], "try-error") ){Rh[[j]]=NA}else{
                       nugget=model[[j]][[type.function]]$nugget
                       Ch=geoR::cov.spatial(H, 
                                            cov.model=model[[j]][[type.function]]$cov.model,
                                            cov.pars=model[[j]][[type.function]]$cov.pars,
                                            kappa=model[[j]][[type.function]]$kappa)+nugget*diag(ncol(H))
                       Rh=Ch/Ch[1,1]
                       Rh[Rh<0]<-0 }
                   }
                   return(Rh)
                 })
  finish_time <- Sys.time()
  stopCluster(cl) # cerrar los cluster
  
  tmp = finish_time - start_time
  print(tmp) #Tiempo de procesamiento
  return(list(Rh=Rh, time=tmp))
}  


## Estimaciones de la matriz de correlación
Rh_np <- function(H, model, type="lambda.opt", type.function="bessel"){
  Rh=list()
  start_time <- Sys.time()
  for(j in 1:length(model)){
    print(paste("Matriz de correlación",j,"de un total de",length(model)))
    if(type=="lambda.opt"){
      C0=model[[j]]$res.estim[[type.function]]$opt$par[1]
      nugget <- line_result2(z=model[[j]]$res.estim[[type.function]]$opt$par, 
                             tj=model[[j]]$res.estim[[type.function]]$tj,
                             uj=model[[j]]$res.estim[[type.function]]$uj, h=0)$gamma
      gamma1=array(NA, dim=dim(H))
      for(i in 1:ncol(H)){
        gamma1[i,] <- line_result2(z=model[[j]]$res.estim[[type.function]]$opt$par, 
                                   tj=model[[j]]$res.estim[[type.function]]$tj,
                                   uj=model[[j]]$res.estim[[type.function]]$uj, h=H[i,])$gamma
        if(gamma1[i,]>C0){gamma1[i,]<-C0}
      }
      Ch <- C0-gamma1+nugget*diag(nrow(H))
      Rh[[j]] <- Ch/Ch[1,1]
    }
    if(type=="cloud"){
      C0=model[[j]][[type.function]]$opt$par[1]
      nugget <- line_result2(z=model[[j]][[type.function]]$opt$par, 
                             tj=model[[j]][[type.function]]$tj,
                             uj=model[[j]][[type.function]]$uj, h=0)$gamma
      gamma1=array(NA, dim=dim(H))
      for(i in 1:ncol(H)){
        gamma1[i,] <- line_result2(z=model[[j]][[type.function]]$opt$par, 
                                   tj=model[[j]][[type.function]]$tj,
                                   uj=model[[j]][[type.function]]$uj, h=H[i,])$gamma
        if(gamma1[i,]>C0){gamma1[i,]<-C0}
      }
      Ch <- C0-gamma1+nugget*diag(nrow(H))
      Rh[[j]] <- Ch/Ch[1,1]
    }
    if(type=="loglik"){
      if(inherits(model[[j]][[type.function]], "try-error") ){Rh[[j]]=NA}else{
        nugget=model[[j]][[type.function]]$nugget
        Ch=geoR::cov.spatial(H, 
                             cov.model=model[[j]][[type.function]]$cov.model,
                             cov.pars=model[[j]][[type.function]]$cov.pars,
                             kappa=model[[j]][[type.function]]$kappa)+nugget*diag(ncol(H))
        Rh[[j]]=Ch/Ch[1,1]
      }
    }
  }
  finish_time <- Sys.time()
  tmp = finish_time - start_time
  print(tmp) #Tiempo de procesamiento
  return(list(Rh=Rh, time=tmp))
}









