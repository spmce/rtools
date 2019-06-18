library(circumplex)

svalSPMCE <- function(data, cvars, evars, cpos=NULL, flip=NULL){
    
    ne <- length(evars)
    nc <- length(cvars)
    
    fflip <- function(x){
        if (any(x > 20 & x < 340)){
            stop("The variables you want to flip are not close enough to 0/360 degrees")
        } else if(all(x <= 20 | x >= 340)) {
            x <- 360 - x
        }
        return(x)
    }
    
    if(is.null(cpos)){
        cpos <- seq(0, 360, by = 360 / nc)[-(nc+1)]
    }
    
    out <- NULL
    
    for(ii in 1:ne){
        temp <- ssm_analyze(.data=data, angles=cpos, scales=cvars, measures=evars[ii], boots=1)
        temp <- as.numeric(temp$results$d_est)
        
        out <- c(out, temp)
    }
    
    if(!is.null(flip)){
        
        if(!is.character(flip) | any((flip %in% evars)==FALSE) | length(flip) > length(evars)){
            stop("'flip' must be a character vector whose elements are a subset of those in 'evars'")
        }
        
        out[evars %in% flip] <- fflip(out[evars %in% flip])
    }
    
    return(out)
}
