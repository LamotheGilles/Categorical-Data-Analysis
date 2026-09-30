 PoissonReg<-function(y,x,epsilon= 1e-8,maxiteration=100)
{
   ## y is a numeric vector of the observed counts
   ## x is the design matrix
   ## n is the number of rows
   n<-length(y)
   ## P is the number of coefficients
   P<-ncol(x)

   ## initial beta (beta0=log(ybar), all other betas are 0)
   ## initial model is the null model 
   coef<-numeric(P)
   coef[1]<-log(mean(y))
   linPred<-x %*% coef
   mu<-as.vector(exp(linPred))

   ## QR-decomposition of Cov(Y)^{1/2}*X
   QR <- qr(diag(sqrt(mu), nrow = length(mu)) %*% x)
   ## inverse of (R'R), which is inverse of Fisher Information
   Cov<- chol2inv(qr.R(QR))
    ## we used QR decompo
    ## alternatively, we can use solve(A)=inverse of A
    ## Cov<-solve(t(x) %*% diag(mu) %*% x) 

   ##
   score<-t(x) %*% matrix((y-mu),ncol=1)


   ## Residual Deviance
   ratio<-y/mu
   log.ratio<-ifelse(ratio==0, 0, log(ratio))
   deviance<-2*sum(y*log.ratio-(y-mu))
   null.dev<-deviance
   df.null<-n-1

   iteration<-0
   test.iteration<-(iteration<maxiteration)
   test.error.large<-TRUE


   while (test.iteration & test.error.large)
   {
    iteration<-iteration+1
    coef<-coef+Cov %*% score

    linPred<-x %*% coef
    mu<-as.vector(exp(linPred))
    ## QR-decomposition of Cov(Y)^{1/2}*X
    QR <- qr(diag(sqrt(mu), nrow = length(mu)) %*% x)
    ## inverse of (R'R), which is inverse of Fisher Information
    Cov<- chol2inv(qr.R(QR))
    ## we used QR decompo
    ## alternatively, we can use solve(A)=inverse of A
    ## Cov<-solve(t(x) %*% diag(mu) %*% x) 
 
   ##
    score<-t(x) %*% matrix((y-mu),ncol=1)

    ratio<-y/mu
    log.ratio<-ifelse(ratio==0, 0, log(ratio))
    dev.old<-deviance
    deviance<-2*sum(y*log.ratio-(y-mu))
    rel.error<-abs(deviance-dev.old)/(abs(deviance)+0.1)
    test.error.large<- (rel.error>=epsilon)

   }


   print(paste("Fisher Scoring Iterations: ",iteration))

   return(list(coefficients=as.vector(coef),vcov<-Cov,deviance=deviance,df=n-P,null.dev=null.dev,df.null=df.null))

}