focalProbpool<-function(x,coords,...) 
{
if(is(x)=="probpool")
  {
  output<-list()
  var.disp.pool<-extract(x@disp.pool,coords)
  var.env.pool<-extract(x@env.pool,coords)
  var.bio.pool<-extract(x@bio.pool,coords)
  output$disp.pools<-var.dist.pool
  output$env.pools<-var.env.pool
  output$bio.pools<-var.bio.pool
  psiDist<-rowSums(var.disp.pool)
  psiEnv <-rowSums(var.env.pool)
  psiBio <-rowSums(var.bio.pool)
  output$psi<-cbind(coordinates(coords),psiDisp,psiEnv,psiBio,psiDisp*psiEnv*psiBio)
  colnames(output$psi)<-c("longitude","latitude","Psi-disp","Psi-env","Psi-bio")
  output
}else{
  stop("x needs to be of class probpool")  
  
  }
  
}

