
setMethod("plot", c("list"),
          function(x, focalunit,...)
            {
              #if(is.vector(x$dist.pool)==TRUE && is.vector(x$env.pool==TRUE && is.vector(x$bio.pool==TRUE)))
              if(length(focalunit)==2){
              par(mfrow=c(3,3))
              plot(x$dist.pool[focalunit[1],focalunit[2]], names(x$dist.pool),  ylab = "Probabilities", xlab= "Species", main = expression(paste(psi, "-dist")))
              barplot(x$dist.pool.pdf[focalunit[1],focalunit[2]],  ylab = "Probabilities", xlab= "Species",, main = expression(paste(psi, "-dist")))
              barplot(rbind(x$dist.pool.cdf[focalunit[1],focalunit[2]], 1-x$dist.pool.cdf[focalunit[1],focalunit[2]]), ylim=c(0, 1),, main = expression(paste(psi, "-dist")))
             
              plot(x$env.pool[focalunit[1],focalunit[2]], names(x$env.pool),  ylab = "Probabilities", xlab= "Species", main = expression(paste(psi, "-env")))
              barplot(x$env.pool.pdf[focalunit[1],focalunit[2]],  ylab = "Probabilities", xlab= "Species", main = expression(paste(psi, "-env")))
              barplot(rbind(x$env.pool.cdf[focalunit[1],focalunit[2]], 1-x$env.pool.cdf[focalunit[1],focalunit[2]]), ylim=c(0, 1), main = expression(paste(psi, "-env")))

              plot(x$bio.pool[focalunit[1],focalunit[2]], names(x$bio.pool),  ylab = "Probabilities", xlab= "Species", main = expression(paste(psi, "-bio")))
              barplot(x$bio.pool.pdf[focalunit[1],focalunit[2]],  ylab = "Probabilities", xlab= "Species", main = expression(paste(psi, "-bio")))
              barplot(rbind(x$env.pool.cdf[focalunit[1],focalunit[2]], 1-x$env.pool.cdf[focalunit[1],focalunit[2]]), ylim=c(0, 1), main = expression(paste(psi, "-bio")))
              }else{
              par(mfrow=c(2,4))  
              raster::plot(x$psi.dist, main=expression(paste(psi, "-dist")))
              raster::plot(x$psi.env,  main=expression(paste(psi, "-env")))
              raster::plot(x$psi.bio,  main=expression(paste(psi, "-bio")))
              raster::plot(x$psi.bio*x$psi.dist*x$psi.env,  main=expression(paste(psi, "-dist+env+bio")))
              raster::plot(x$psi.env*x$psi.dist,  main=expression(paste(psi, "-dist+env")))
              raster::plot(x$psi.bio*x$psi.dist,  main=expression(paste(psi, "-dist+bio")))
              raster::plot(x$psi.bio*x$psi.env,  main=expression(paste(psi, "-env+bio")))
              
              }
          }
)




