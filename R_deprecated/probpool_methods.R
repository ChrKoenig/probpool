
# set standard method for plotting a class probpool

setMethod("plot", c("probpool"),
          function(x, focalunit=FALSE,...)
            {
              if(length(focalunit)==2)
                {
                par(mfrow=c(4,x@pool.count))
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
                if(is.null(x@pools$comb.pool)==FALSE)
                  {
                    plotPoolProbs(x,"comb.pool",focalunit)
                    plotPoolPDF(x,"comb.pool",focalunit)
                    plotPoolCDF(x,"comb.pool",focalunit)
                  }
              }else{
              par(mfrow=c(2,x@pool.count))    
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
                if(is.null(x@PSI$comb.pool)==FALSE)  
                  {
                    plotRasterPool(x,"comb.pool")
                  }
              }
          }
)

