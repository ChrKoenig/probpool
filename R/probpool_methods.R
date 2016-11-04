
setMethod("plot", c("probpool"),
          function(x, focalunit=FALSE,...)
          {
            if(length(focalunit)==2){
              par(mfrow=c(4,3))
              focal<-cbind(focalunit[1],focalunit[2])
              loc<-extract(x@pools$disp.pool,focal, cellnumbers=TRUE)
              
              barplot(rev(sort(x@pools$disp.pool[loc[1]][1,])), ylim=c(0,1), ylab = "Probabilities", xlab= "Species", main = expression(paste(psi, "-dist")))
              barplot(dpoibin(seq(1,length(x@pools$disp.pool[focalunit[1],focalunit[2]][1,]),1), x@pools$disp.pool[focalunit[1],focalunit[2]][1,]), ylab = "Probabilities", xlab= "Species", main = expression(paste(psi, "-dist")))
              barplot(ppoibin(seq(1,length(x@pools$disp.pool[focalunit[1],focalunit[2]][1,]),1), x@pools$disp.pool[focalunit[1],focalunit[2]][1,]), ylim=c(0,1), ylab = "Probabilities", xlab= "Species", main = expression(paste(psi, "-dist")))
              
               
              barplot(rev(sort(x@pools$env.pool[loc[1]][1,])), ylim=c(0,1), ylab = "Probabilities", xlab= "Species", main = expression(paste(psi, "-env")))
              barplot(dpoibin(seq(1,length(x@pools$env.pool[loc[1]][1,]),1), x@pools$env.pool[focalunit[1],focalunit[2]][1,]), ylab = "Probabilities", xlab= "Species", main = expression(paste(psi, "-dist")))
              barplot(ppoibin(seq(1,length(x@pools$env.pool[loc[1]][1,]),1), x@pools$env.pool[focalunit[1],focalunit[2]][1,]), ylim=c(0,1), ylab = "Probabilities", xlab= "Species", main = expression(paste(psi, "-dist")))
              
                 
              barplot(rev(sort(x@pools$bio.pool[loc[1]][1,])), ylim=c(0,1), ylab = "Probabilities", xlab= "Species", main = expression(paste(psi, "-bio")))
              barplot(dpoibin(seq(1,length(x@pools$bio.pool[loc[1]][1,]),1), x@pools$bio.pool[focalunit[1],focalunit[2]][1,]), ylab = "Probabilities", xlab= "Species", main = expression(paste(psi, "-dist")))
              barplot(ppoibin(seq(1,length(x@pools$bio.pool[loc[1]][1,]),1), x@pools$bio.pool[focalunit[1],focalunit[2]][1,]), ylim=c(0,1), ylab = "Probabilities", xlab= "Species", main = expression(paste(psi, "-dist")))
              
             
              barplot(rev(sort(x@pools$bio.pool[loc[1]][1,]*x@pools$disp.pool[focalunit[1],focalunit[2]][1,]*x@pools$env.pool[focalunit[1],focalunit[2]][1,])), ylim=c(0,1), ylab = "Probabilities", xlab= "Species", main = expression(paste(psi, "disp*env*bio")))
              barplot(dpoibin(seq(1,length(x@pools$bio.pool[loc[1]][1,]),1), x@pools$bio.pool[loc[1]][1,]*x@pools$disp.pool[loc[1]][1,]*x@pools$env.pool[loc[1]][1,]), ylab = "Probabilities", xlab= "Species", main = expression(paste(psi, "-dist")))
              barplot(ppoibin(seq(1,length(x@pools$bio.pool[loc[1]][1,]),1), x@pools$bio.pool[loc[1]][1,]*x@pools$disp.pool[loc[1]][1,]*x@pools$env.pool[loc[1]][1,]), ylim=c(0,1), ylab = "Probabilities", xlab= "Species", main = expression(paste(psi, "-dist")))
              
              fullpool.rst<-x@pools$disp.pool*x@pools$env.pool*x@pools$bio.pool
              fullpool    <-raster::calc(fullpool.rst,fun=sum)
              
              }else{
              par(mfrow=c(2,4))  
              raster::plot(x@PSI$PSI.disp, main=expression(paste(psi, "-disp")))
              raster::plot(x@PSI$PSI.env,  main=expression(paste(psi, "-env")))
              raster::plot(x@PSI$PSI.bio,  main=expression(paste(psi, "-bio")))
              raster::plot(fullpool,  main=expression(paste(psi, "-disp+env+bio")))
              
              
            }
          }
)

# Function to plot pool for a single focal unit defined by focalunit
plotPoolProbs<-function(probpool,pool,focalunit)
{
  focal<-cbind(focalunit[1],focalunit[2])
  loc<-extract(probpool@pools[[pool]],focal, cellnumbers=TRUE)
  barplot(rev(sort(probpool@pools[[pool]][loc[1]][1,])), ylim=c(0,1), ylab = "Probabilities", xlab= "Species", main = pool)
}

plotPoolPDF<-function(probpool,pool,focalunit)
{
  focal<-cbind(focalunit[1],focalunit[2])
  loc<-extract(probpool@pools[[pool]],focal, cellnumbers=TRUE)
  barplot(dpoibin(seq(1,length(probpool@pools[[pool]][focalunit[1],focalunit[2]][1,]),1), probpool@pools[[pool]][focalunit[1],focalunit[2]][1,]), ylab = "Probabilities", xlab= "Species", main = expression(paste(psi, "-dist")))
}


# Function to plot pool as a raster
plotRasterPool<-function(probpool,pool)
{
  raster::plot(probpool@PSI[[pool]], main=pool)
}


