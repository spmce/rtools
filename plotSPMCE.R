library(ggplot2)

plotSPMCE <- function(res, clabels=NULL, elabels=NULL, print.pdfs=TRUE, keep.plots=FALSE){
    
    # PLOT 1: CIRCUMPLEX PLOT WITH COVARIATE PROJECTIONS AND CONFIDENCE ELLIPSES
    # HELPERS
    drawcirc <- function(center = c(0,0), diameter=2, npoints=500){
        r = diameter / 2
        tt <- seq(0,2*pi,length.out = npoints)
        xx <- center[1] + r * cos(tt)
        yy <- center[2] + r * sin(tt)
        return(data.frame(x = xx, y = yy))
    }
    
    cr1 <- drawcirc(diameter = 2, npoints=10000)
    
    pos <- data.frame(
        cos = cos(res$theta$est),
        sin = sin(res$theta$est)
    )
    
    # SET LABELS
    if(is.null(clabels)){
        clabels <- rownames(res$theta)
    }
    
    # BUILD BASIC CIRCUMPLEX
    pc <- ggplot()
    
    pc <- pc + 
        geom_segment(aes(x=1,y=0,xend=-1,yend=0), colour="darkgrey", size=.5)+
        geom_segment(aes(x=0,y=-1,xend=0,yend=1), colour="darkgrey", size=.5)+
        geom_point(aes(x=0,y=0), size=.75, colour="darkgrey")+
        geom_path(aes(x, y), data = cr1, lwd=1.05)+
        #geom_point(aes(x=cos,y=sin), data=pos, shape=15, size=3)+
        scale_y_continuous(limits=c(-1.25, 1.25))+
        scale_x_continuous(limits=c(-1.25, 1.25))+
        geom_text(aes(x=cos*(1.15), y=sin*(1.15)), data=pos, label=clabels, size=3, parse=TRUE)+
        geom_spoke(aes(x=cos, y=sin, angle=res$theta$est, radius=.075), data=pos, size=1)
    
    # CHECK FOR COVARIATES
    ncov <- length(res$covs)
    
    # SET LABELS
    if(is.null(elabels)){
        elabels <- paste0("EV",1:ncov)
    }

    # COVARIATES AND CONFIDENCE ELLIPSES
    if(ncov > 0){
        
        phi <- seq(0, 360, length.out=500)
        
        points <- list()
        ellipses <- list()
        
        for(ii in 1:ncov){
            
            temp <- data.frame(
                pc = res$covs[[ii]][paste0("PHI_C", LETTERS[ii]),"est"],
                ps = res$covs[[ii]][paste0("PHI_S", LETTERS[ii]),"est"],
                a = res$covci[[ii]]["a_q", "est"],
                b = res$covci[[ii]]["b_q", "est"],
                w = res$covci[[ii]]["w_q", "est"]
            )
                
            ellipses[[paste0("EV",ii)]] <- data.frame(
                cos = with(cbind(temp,phi), pc+a*cos(w)*cos(phi*pi/180)-b*sin(w)*sin(phi*pi/180)),
                sin = with(cbind(temp,phi), ps+a*sin(w)*cos(phi*pi/180)+b*cos(w)*sin(phi*pi/180))
            )
            
            points[[paste0("EV",ii)]] <- data.frame(x=temp$pc, y=temp$ps)
            
            pc <- pc + 
                geom_point(data=points[[ii]], aes(x,y), size=1.05)+
                geom_text(data=points[[ii]], aes(x, y+.05), size=2, label=elabels[ii], parse=TRUE)+
                geom_path(data=ellipses[[ii]], aes(x=cos, y=sin), linetype="dashed", size=.33)
        }
        
            
    
        
    }
    
    # FINISHING TOUCH
    pc <- pc+
        theme(panel.background = element_blank())+
        theme(panel.grid.major = element_blank())+
        theme(panel.grid.minor = element_blank())+
        theme(axis.line = element_blank())+
        theme(axis.ticks = element_blank())+
        theme(axis.text.x = element_blank())+
        theme(axis.text.y = element_blank())+
        theme(axis.title.x = element_blank())+
        theme(axis.title.y = element_blank())+
        theme(plot.title = element_text(size=20, hjust=.5))+
        theme(legend.position="none")
    
    # PLOT 2: CORRELATION FUNCTION
    dcf <- data.frame(dphi=seq(0,360,by=.5))
    
    corf <- function(b, dphi){
        
        nbm <- length(b[-1])
        
        rho <- rep(b[1]^2, length(dphi))
        
        for(ii in 1:nbm){
            rho <- rho + b[ii+1]^2 * cos(ii*dphi) 
        }
        
        return(rho)
    }
    
    dcf$rho <- corf(b=res$beta$est, dphi=(dcf$dphi*pi)/180)
    
    pcf <- ggplot(dcf)+
        geom_hline(yintercept=0, colour="darkgrey", size=.5)+
        geom_line(aes(dphi,rho))+
        scale_y_continuous(limits=c(-1.1,1.1))+
        scale_x_continuous(limits=c(0,360),
                           breaks=seq(0,360,45), 
                           labels=parse(text=paste0(seq(0,360,45), "*degree")))+
        labs(y=expression(Common~Score~Correlations~rho[paste(c[j], ",", c[k])]),
             x=expression(Separating~Angle~(theta[j]-theta[k])))+
        theme_classic()+
        theme(plot.title = element_text(size=20, hjust=.5))+
        theme(legend.position="none")+
        theme(axis.ticks=element_line(color="black"))+
        theme(axis.text=element_text(color="black"))
    
    # PLOT(S) 3: PROFILE PLOTS
    if(ncov > 0){
        
        pcp <- list()
        
        dfcp <- list()
        dfcpu <- list()
        
        tu <- res$unq$est
        tb <- as.numeric(res$beta$est)
        tt <- res$theta$est
        
        nbm <- length(tb[-1])
        
        dcp <- data.frame(phi=seq(0,360,by=.5))
        
        for(ii in 1:ncov){
            
            tur <- res$r_unq[[ii]]$est
            tg0 <- res$covs[[ii]][paste0("G0",LETTERS[ii]), "est"]
            tg1 <- res$covs[[ii]][paste0("G1",LETTERS[ii]), "est"]
            td <- res$covs[[ii]][paste0("D",LETTERS[ii]), "est"]

            # FUNCTION
            rho <- tg0 * rep(tb[1]^2, length(dcp$phi))
            
            t1 <- 0
            
            for(jj in 1:nbm){
                t1 <- t1 + tb[jj+1]^2 * cos(jj*((dcp$phi*pi)/180 - td))
            }
            
            rho <- rho + tg1 * t1
            
            # WITH ERROR
            rho_u <- tg0 * rep(tb[1]^2, length(tt))
            
            t1 <- 0
            
            for(jj in 1:nbm){
                t1 <-  t1 + tb[jj+1]^2 * cos(jj*(tt - td)) 
            }
            
            rho_u <- (rho_u + tg1 * t1) + tur*tu
            rho_u <- c(rho_u, rho_u[1])
            
            lims <- round((range(c(rho_u, rho)) * 1.25) + c(-.1,.1), 1)
            
            # DF FOR GGPLOT
            dfcp[[ii]] <- data.frame(phi=dcp$phi, rho=rho)
            dfcpu[[ii]] <- data.frame(phi=c((tt*180)/pi, 360), rho=rho_u)
            
        }
        
        if(ncov==1){
            lims <- round(range(c(rho_u, rho) * 1.25) + c(-.1,.1), 1)
        } else if(ncov>1) {
            lims <- round(range(c(do.call("rbind", dfcp)$rho, 
                                  do.call("rbind", dfcpu)$rho))* 1.25 + c(-.1,.1), 1)
        }
        
        for(ii in 1:ncov){
            
            pp <- ggplot()+
                geom_hline(yintercept=0, colour="darkgrey", size=.5)+
                geom_line(data=dfcp[[ii]],aes(x=phi, y=rho))+
                geom_point(data=dfcpu[[ii]], aes(x=phi, y=rho), pch=1, size=2)+
                scale_y_continuous(limits=lims, breaks=round(seq(lims[1], lims[2], .1),1))+
                scale_x_continuous(limits=c(0,360),
                                   breaks=seq(0,360,45), 
                                   labels=parse(text=paste0(seq(0,360,45), "*degree")))+
                labs(y="External Correlation Function (Solid Line)\nand External Covariances (Dots)",
                     x=expression(Angular~Position~(theta[j])))+
                ggtitle(substitute(paste("Covariate: ", elab), list(elab=elabels[[ii]])))+
                theme_classic()+
                theme(plot.title = element_text(size=14, hjust=.5))+
                theme(axis.ticks=element_line(color="black"))+
                theme(axis.text=element_text(color="black"))
            
            pcp[[ii]] <- pp
            
        }
        
    }
    
    # OUTPUT
    out <- list(
        circplot = pc,
        corfun = pcf
    )
    
    if(ncov > 0){
        out <- c(out,
                 profplot = pcp)
    }
    
    # CREATE PDFS
    if(print.pdfs==TRUE){
        filenames <- paste0(gsub("\\.out","",res$file), "_", names(out), ".pdf")
    
        for(ii in 1:length(filenames)){
        
            pdf(file=filenames[ii], height=6, width=6)
            print(out[[ii]])
            dev.off()
        }
    }
    
    # RETURN PLOTS?
    if(keep.plots==TRUE){
        return(out)
        
        for(ii in 1:length(out)) print(out[[ii]])
        
    }
}
