
setMethod("plot", c("probpool"),
          function(x, focalunit=FALSE,...)
          {
            if(length(focalunit)==2)
              {
              par(mfrow=c(4,3))
              
              if(is.null(x@pools$env.pool)==FALSE)
              {
                plotPoolProbs(x,"env.pool",focalunit)
                plotPoolPDF(x,"env.pool",focalunit)
                plotPoolCDF(x,"env.pool",focalunit)
              }
              
              if(is.null(x@pools$disp.pool)==FALSE)
              {
                plotPoolProbs(x,"disp.pool",focalunit)
                plotPoolPDF(x,"disp.pool",focalunit)
                plotPoolCDF(x,"disp.pool",focalunit)
              }
              
              if(is.null(x@pools$bio.pool)==FALSE)
              {
                plotPoolProbs(x,"bio.pool",focalunit)
                plotPoolPDF(x,"bio.pool",focalunit)
                plotPoolCDF(x,"bio.pool",focalunit)
              }
              
              #fullpool.rst<-x@pools$disp.pool*x@pools$env.pool*x@pools$bio.pool
              #fullpool    <-raster::calc(fullpool.rst,fun=sum)
              
              }else{
                
              par(mfrow=c(2,4))    
                
                if(is.null(x@PSI$bio.pool)==FALSE)  
                {
                  plotRasterPool(x,"bio.pool")
                }
                
                if(is.null(x@PSI$env.pool)==FALSE)  
                {
                  plotRasterPool(x,"env.pool")
                }
                
                if(is.null(x@PSI$bio.pool)==FALSE)  
                {
                  plotRasterPool(x,"disp.pool")
                }

              #raster::plot(fullpool,  main=expression(paste(psi, "-disp+env+bio")))

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

# Function to plot pool probability density function
plotPoolPDF<-function(probpool,pool,focalunit)
{
  focal<-cbind(focalunit[1],focalunit[2])
  loc<-extract(probpool@pools[[pool]],focal, cellnumbers=TRUE)
  barplot(dpoibin(seq(1,length(probpool@pools[[pool]][loc[1]][1,]),1), probpool@pools[[pool]][loc[1]][1,]), names.arg=1:length(probpool@pools[[pool]][loc[1]][1,]), ylab = "Probabilities", xlab= "Species", main = pool)
}

# Function to plot pool cumulative density function
plotPoolCDF<-function(probpool,pool,focalunit)
{
  focal<-cbind(focalunit[1],focalunit[2])
  loc<-extract(probpool@pools[[pool]],focal, cellnumbers=TRUE)
  barplot(ppoibin(seq(1,length(probpool@pools[[pool]][loc[1]][1,]),1), probpool@pools[[pool]][loc[1]][1,]), ylim=c(0,1), names.arg=1:length(probpool@pools[[pool]][loc[1]][1,]), ylab = "Probabilities", xlab= "Species", main = pool)
}

# Function to plot pool as a raster
plotRasterPool<-function(probpool,pool)
{
  raster::plot(probpool@PSI[[pool]], main=pool)
}


