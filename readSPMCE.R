readSPMCE <- function(path, conf=.05) {
    
    if(sum(grepl("NO CONVERGENCE", readLines(path))) > 0){
        stop("\nThe Mplus algorithm did not converge. 
             Try adjusting the starting values pertaining to the positions.")
    }

    out <- readModels(path)
    
    if(length(out$errors) > 0){
        stop("\nIt appears as if Mplus did not find an appropriate solution. 
             Try adjusting the starting values pertaining to the positions.")
    }
    
    # FIT INFORMATION
    fit <- readLines(path)
    
    start <- which(grepl(fit, pattern="MODEL FIT INFORMATION"))
    end <- which(grepl(fit, pattern="MODEL RESULTS")) - 1
    
    fit <- fit[start:end]
    
    class(fit) <- "mpfit"
    
    # PARAMETER ESTIMATES
    pars <- out$parameters$unstandardized
    
    pars <- pars[pars$pval!=999,]
    
    # NUMBER OF COVARIATES
    ncov <- length(unique(pars[grepl("^L_EV",pars$param),"param"]))
    
    # OUTPUT
    res <- list()
    
    # FILENAME
    res$file <- path
    
    # FIT INFORMATION
    res$fit <- fit
    
    # ZETA
    res$zeta <- pars[grep(pars$paramHeader, pattern="^L_C.*BY$"),c("est","se","pval")]
    
    # NU
    res$unq <- pars[grep(pars$paramHeader, pattern="^U_C.*BY$"),c("est","se","pval")]
    
    # THETA
    res$theta <- rbind(c(0,0,999),
                   pars[grepl(pars$paramHeader, pattern="^New") & 
             grepl(pars$param, pattern="^T"), 
         c("est","se","pval")])
    
    res$theta[,"est (deg)"] <- (res$theta$est * 180) / pi
    res$theta[,"se (deg)"] <- (res$theta$se * 180) / pi
    
    # COMMUNALITIES
    res$h <- (1+res$unq[,1,drop=FALSE]^2)^-(1/2)
    
        rownames(res$zeta) <-  
        rownames(res$unq) <- 
        rownames(res$theta) <- 
        rownames(res$h) <- paste0("CV", 1:nrow(res$zeta))

    # BETAS
    ncv <- nrow(res$zeta)
    ncos <- nrow(pars[grep(pars$paramHeader, pattern="^COS.*BY$"),c("est","se","pval")])
    fcomp <- ncos/ncv
    
    res$beta <- pars[grep(pars$paramHeader, pattern="^LEV.*BY$")[1],c("est","se","pval")]
    
    for(ii in 1:fcomp){
        res$beta <- rbind(res$beta,
                          pars[grep(pars$paramHeader, 
                                    pattern=paste0("^COS", ii, ".*BY$"))[1],
                               c("est","se","pval")])
    }
    
    rownames(res$beta) <- c("beta_0",paste0("beta_",1:fcomp))
    
    # COVARIATE INFORMATION
    if(ncov > 0){
        
        chiv <- qchisq(df=2, p=(1-conf))

        tech3 <- out$tech3$paramCov
        
        m <- list()
        
        idx <- ((nrow(tech3)-2*ncov)+1):nrow(tech3)
        idx <- split(idx, ceiling(seq_along(idx)/2))
        
        # PARAMETER ESTIMATES
        res$covs <- list()
        
        for(ii in 1:ncov){
            res$covs[[paste0("Covariate: EV", ii)]] <- 
                pars[grepl(pars$param, pattern=paste0(LETTERS[ii], "$")), 
                     c("est","se","pval")]
            
            rownames(res$covs[[paste0("Covariate: EV", ii)]]) <- 
                pars[grepl(pars$param, pattern=paste0(LETTERS[ii], "$")), 
                     c("param")]
            
            res$covs[[ii]] <- rbind(res$covs[[ii]],
                                    (cbind(res$covs[[ii]][paste0("D",LETTERS[ii]),1:2], pval=NA)*180)/pi
            )
            
            rownames(res$covs[[ii]])[nrow(res$covs[[ii]])] <- paste0("D", LETTERS[ii], " (deg)")
            }
        
        # INTERVAL ESTIMATES
        for(ii in 1:ncov){
        
            m[[ii]] <- tech3[idx[[ii]], idx[[ii]]]
            m[[ii]][1,2] = m[[ii]][2,1]
        
            meig <- sqrt(eigen(m[[ii]] * chiv)$values)
        
            phi_c <- res$covs[[ii]][paste0("PHI_C",LETTERS[ii]),"est"]
            phi_s <- res$covs[[ii]][paste0("PHI_S",LETTERS[ii]),"est"]
        
            a <- meig[1]
            b <- meig[2]
        
            vc <- m[[ii]][1,1]
            vs <- m[[ii]][2,2]
            ccs <- m[[ii]][2,1]
            rcs <- ccs/(sqrt(vc*vs))
        
            if(ccs > 0 & vc > vs){
                w <- .5 * (atan(2*ccs / (vc-vs)))
            }
        
            if(vc < vs){
                w <- .5 * (atan(2*ccs / (vc-vs))+pi)
            }
        
            if(ccs < 0 & vc > vs){
                w <- .5 * (atan(2*ccs / (vc-vs))+2*pi)
            }
            
            area <- pi * chiv * sqrt(vs) * sqrt(vc) * (sqrt(1-rcs^2))
            
            phi <- seq(0, 360, length.out=1000)
            
            coords <- data.frame(c=phi_c+a*cos(w)*cos(phi)-b*sin(w)*sin(phi),
                                 s=phi_s+a*sin(w)*cos(phi)+b*cos(w)*sin(phi)
            )

            res$covci[[paste0("Covariate: EV", ii)]] <-
                data.frame(est=c(a, b, w, area))
            
            rownames(res$covci[[paste0("Covariate: EV", ii)]]) <- 
                c("a_q","b_q","w_q", "area")
            
            originIncluded <- nrow(unique(sign(coords))) == 4
            
            # INTERVAL ESTIMATES SENSITIVITY
            if(!originIncluded){
                
                vq <- (0-(phi_c))*cos(w)+(0-phi_s)*sin(w)
                wq <- -(0-(phi_c))*sin(w)+(0-phi_s)*cos(w)
            
                fun <- function(t){
                    -wq*b+(2*vq*a-2*(a^2-b^2))*t+(2*vq*a+2*(a^2-b^2))*t^3+wq*b*t^4
                }
                
                h1 <- a*b*vq*wq*(a^2-b^2)
                h2 <- (((a^2-b^2)^2-vq^2*a^2-wq^2*b^2)/3)^3
                
                Math.cbrt <- function(x){
                    sign(x) * abs(x)^(1/3)
                }
                
                h3 <- Math.cbrt(sqrt(h1^2-h2)-h1)
                
                A = (2*h3)/(wq*b)+(2*((a^2-b^2)^2-vq^2*a^2-wq^2*b^2)) / 
                    (3*wq*b*h3)

                B = ((vq*a+(a^2-b^2))/(wq*b))^2
                
                C = -8*(((vq*a+(a^2-b^2))/(wq*b))^3+2*((vq*a-(a^2-b^2))/(wq*b)))
                
                D = -(vq*a+a^2-b^2)/(2*wq*b)
                
                roots <- suppressWarnings(c(
                (-1/2*sqrt(A+B))+(-1/2*sqrt(2*B-A-C/(4*sqrt(A+B)))) + D,
                
                (-1/2*sqrt(A+B))+(1/2*sqrt(2*B-A-C/(4*sqrt(A+B)))) + D,
                
                (1/2*sqrt(A+B))+(-1/2*sqrt(2*B-A+C/(4*sqrt(A+B)))) + D,

                (1/2*sqrt(A+B))+(1/2*sqrt(2*B-A+C/(4*sqrt(A+B)))) + D
                ))
                
                roots <- sort(as.numeric(na.omit(roots)))

                r1 <- roots[1]
                r2 <- roots[2]
        
                cp1 <- (1-r1^2)/(1+r1^2)
                sp1 <- (2*r1) / (1+r1^2)
                cp2 <- (1-r2^2)/(1+r2^2)
                sp2 <- (2*r2) / (1+r2^2)
        
                dists <- c(sqrt((vq-a*cp1)^2 + (wq-b*sp1)^2),
                           sqrt((vq-a*cp2)^2 + (wq-b*sp2)^2))
        
                dmin <- min(dists)
                dmax <- max(dists)

                # INTERVAL ESTIMATES DISPLACEMENT
                m1 <- ((phi_c*phi_s-(ccs*chiv)) / (phi_c^2-(vc*chiv))) +
                    ((chiv*(vs*(phi_c^2-(1-rcs^2)*vc*chiv)+phi_s^2*vc-2*ccs*phi_c*phi_s))^(1/2)) / 
                    (phi_c^2-vc*chiv)

                m2 <- ((phi_c*phi_s-(ccs*chiv)) / (phi_c^2-(vc*chiv))) -
                    ((chiv*(vs*(phi_c^2-(1-rcs^2)*vc*chiv)+phi_s^2*vc-2*ccs*phi_c*phi_s))^(1/2)) / 
                    (phi_c^2-vc*chiv)
                
                c1 <- (vs*phi_c+vc*m1*phi_s-ccs*(phi_s+m1*phi_c)) /
                    (vs+m1^2*vc-2*m1*ccs)
                
                c2 <- (vs*phi_c+vc*m2*phi_s-ccs*(phi_s+m2*phi_c)) /
                    (vs+m2^2*vc-2*m2*ccs)

                s1 <- c1 * m1
                s2 <- c2 * m2
                
                if(c1 > 0 & s1 > 0){
                    d1 <- atan(s1/c1)
                }
                
                if(c1 < 0){
                    d1 <- atan(s1/c1) + pi
                }
                
                if(c1 > 0 & s1 < 0){
                    d1 <- atan(s1/c1) + 2*pi
                }

                if(c2 > 0 & s2 > 0){
                    d2 <- atan(s2/c2)
                }
                
                if(c2 < 0){
                    d2 <- atan(s2/c2) + pi
                }
                
                if(c2 > 0 & s2 < 0){
                    d2 <- atan(s2/c2) + 2*pi
                }
                
                res$covci[[paste0("Covariate: EV", ii)]] <- 
                    rbind(res$covci[[paste0("Covariate: EV", ii)]],
                    data.frame(est=c(dmin, dmax, (d1*180)/pi, (d2*180)/pi))
                    )
                
                rownames(res$covci[[paste0("Covariate: EV", ii)]]) <- 
                    c(rownames(res$covci[[paste0("Covariate: EV", ii)]])[1:4], 
                      c("gamma_1q", "gamma_2q", "d1 (deg)","d2 (deg)"))
            }
            
            # UNIQUE CORRELATIONS
            res$r_unq <- list()
            
            for(ii in 1:ncov){
                res$r_unq[[paste0("Covariate: EV", ii)]] <- 
                    pars[grepl(pars$paramHeader, pattern="^U_.*WITH$") & 
                             grepl(pars$param, pattern=paste0("L_EV", ii)), 
                         c("est","se","pval")]
                
                rownames(res$r_unq[[paste0("Covariate: EV", ii)]]) <- paste0("CV", 1:nrow(res$zeta))
            }
            
            # RSQUARE MEASURES
            res$RSq <- list()
            
            tu <- res$unq$est
            tb <- as.numeric(res$beta$est)
            tt <- res$theta$est
            nbm <- length(tb[-1])
            
            for(ii in 1:ncov){
                
                tur <- res$r_unq[[ii]]$est
                tg0 <- res$covs[[ii]][paste0("G0",LETTERS[ii]), "est"]
                tg1 <- res$covs[[ii]][paste0("G1",LETTERS[ii]), "est"]
                td <- res$covs[[ii]][paste0("D",LETTERS[ii]), "est"]
                
                rho_u <- tg0 * rep(tb[1]^2, length(tt))
                
                t1 <- 0
                
                for(jj in 1:nbm){
                    t1 <-  t1 + tb[jj+1]^2 * cos(jj*(tt - td)) 
                }
                
                rho <- (rho_u + tg1 * t1)
                rho_u <- (rho_u + tg1 * t1) + tur*tu
                
                RPq <- cor(rho, rho_u)^2
                
                gt <- res$covs[[ii]]
                
                g0 <- gt[grepl("^G0", rownames(gt)), "est"]
                g1 <- gt[grepl("^G1", rownames(gt)), "est"]
                
                RFq <- g0^2 * res$beta$est[1]^2 + g1^2 * (1-g0^2) 
                RDq <- sum(res$r_unq[[ii]]$est^2) 
                RTq <- RFq + RDq
                
                res$RSq[[paste0("Covariate: EV", ii)]] <- data.frame(
                    est=c(RFq, RDq, RTq, RPq)
                )
                
                rownames(res$RSq[[paste0("Covariate: EV", ii)]]) <-
                    c("R2_Fq", "R2_Dq", "R2_Tq", "r2_Pq")
            }
        }
    }
    
    # WALD TESTS
    if(any(grepl("Wald", names(out$summaries)))){
        res$wald <- out$summaries[grepl("Wald", names(out$summaries))]
    }
    
    return(res)

}

# PRINT METHOD FOR FIT INFORMATION
print.mpfit <- function(x,...){cat(x, sep="\n")}





