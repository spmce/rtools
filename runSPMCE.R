runSPMCE <- function(
  data, 
  file,
  cvars,
  evars=NULL,
  fcomp=1,
  start.cpos=NULL,
  start.epos=NULL,
  eq.scaling=FALSE,
  eq.unique=FALSE,
  pos.fixed=FALSE,
  wald.test=c("none", "zero_sens", "coordinate", "sensitivity", "displacement"),
  wald.vars=evars,
  run=TRUE,
  miss=-99){
    
  
    # INPUT CHECKS AND HELPERS
    if(!is.data.frame(data)) stop("Input 'data' must be an object a dataframe.")
  
    if(fcomp > 1){
        betas <- round(1/(exp(1:fcomp)),1)+.1
    } else {
        betas <- .6
    }
    
    nc <- length(cvars)
    ne <- length(evars) 
    
    if(nc/2 < fcomp) stop("Number of Fourier Components must be less than the number of circumplex indicators divided by 2!")
    
    # STARTING VALUES FOR COVARIATE INTERCORRELATIONS
    if(ne > 1){
        
        svalCovariateCor <- cor(data[, evars], use="pairwise.complete.obs")
        svalCovariateCor <- svalCovariateCor[upper.tri(svalCovariateCor,diag=FALSE)]
        svalCovariateCor <- round(svalCovariateCor,2)
        
        svalCovariateCor <- data.frame(sval=svalCovariateCor, v1=NA, v2=NA)
        
        svalCovariateCor [,"v1"] <- rep(1:(ne-1), times=((ne-1):1)) 
        svalCovariateCor [,"v2"] <- matrix(1:ne, ncol=ne, nrow=ne)[lower.tri(matrix(1:ne, ncol=ne, nrow=ne))]
    }
    
    # SET MISSING DATA
    if(any(is.na(data))){
        message(paste0("Note: NAs were recoded to ", miss))
        }
    
    data[is.na(data)] <- miss
    
    # SET NAMES
    if(ne > 0){
        x <- data[,(c(cvars, evars))]
        colnames(x) <- c(paste0("CV",1:nc), paste0("EV",1:ne))
    } else {
        x <- data[,(cvars)]
        colnames(x) <- paste0("CV",1:nc)
    }
    
    vc <- colnames(x)[1:nc]
    ve <- colnames(x)[(nc+1):(ne+nc)]
    
    # STARTING VALUES FOR CIRCUMPLEX INDICATOR DISCPLACEMENTS
    if(is.null(start.cpos)){
        start.cpos <- seq(0, 360, by = 360 / nc)[2:nc]
    } else {
        start.cpos <- start.cpos[2:nc]
    }
    
    # STARTING VALUES FOR COVARIATE DISPLACEMENTS
    if(is.null(start.epos)){
        start.epos <- sample(seq(0,2*pi,.1), size=ne, replace=TRUE)
    } else if(length(start.epos)!=ne) {
      stop("Argument 'start.epos' must have the same length as argument 'evars'")
    } else {
        nna <- sum(is.na(start.epos))
        start.epos <- (start.epos*pi) / 180
        start.epos[is.na(start.epos)] <- sample(seq(0,2*pi,.1), size=nna, replace=TRUE)
    }
    
    # CREATING INPUT FILE TITLE
    title <- c("\t\t\tCircumplex Indicators:", 
                   paste0("\t\t\tCV", 1:nc, ": ", cvars[1:nc]))
    if(ne > 0){
        title <- c(title, "", "\t\t\tCovariates:", 
                   paste0("\t\t\tEV", 1:ne, ": ", evars[1:ne]))
    }
        
    title <- c(title, "", paste0("\t\t\tNumber of Fourier Components: ", fcomp, ";"))
    
    # WRITE DATA
    write.table(x, file=paste0(file,".dat"), 
                col.names=FALSE, 
                row.names=FALSE, 
                sep="\t", 
                dec=".", 
                quote=FALSE)
    
    # CREATE INPUT FILE
    file.create(paste0(file,".inp"))
    
    # CREATE SYNTAX FOR INPUT FILE
    inpM <- c("TITLE:\t\tSPMC(E) SETUP","", title, "",
                   paste0("DATA: \t\tFILE IS ", file, ".dat", ";"), "")
    
    inpM <- c(inpM, "VARIABLE: \tNAMES ARE",
    
    paste0("\t\t\t", c(colnames(x), ";")), "",
    
    "\t\t\tUSEVARIABLES ARE",
    
    paste0("\t\t\t", c(colnames(x), ";")), "",
    
    paste0("\t\t\tMISSING ARE ALL (", miss, ");"), "",
    
    "ANALYSIS:\t TYPE IS GENERAL;",
    
    "\t\t\tITERATIONS ARE 10000;",
    
    "\t\t\tPROCESSORS ARE 4;", 
    
    "",
    
    "MODEL:\t\t!MEASUREMENT PART FOR THE CIRCUMPLEX",
    
    "\t\t\t!SCALING CONSTANTS", "",
    
    paste0(paste0("\t\t\t", "L_", vc), " BY ", vc , "*1 (L", 1:nc, ");"), "",
    
    "\t\t\t!UNIQUE FACTOR LOADINGS:", "",
    
    paste0(paste0("\t\t\t","U_", vc), " BY ", paste0("L_", vc), "*.5 (U", 1:nc, ");"), "",
    
    "\t\t\t!COMMON FACTORS:", "",
    
    paste0("\t\t\tLEV BY ", paste0("L_", vc[1]), "-", paste0("L_", vc[nc]), "*.4 (b0);"),
    
    "")
    
    for(ii in 1:fcomp){
        inpM <- c(inpM, 
                  paste0("\t\t\tCOS", ii, " BY L_", vc[1], "*", betas[ii], " (b", ii, ");"),
                  paste0("\t\t\tCOS",ii, " BY L_", vc[2:nc], "*", 
                         round(cos(start.cpos * (pi/180)),3),
                         " (C", 2:nc, ii, ");"), 
                  "", 
                  paste0("\t\t\tSIN", ii, " BY L_", vc[1], "@0;"),
                  paste0("\t\t\tSIN",ii, " BY L_", vc[2:nc], "*", 
                         round(sin(start.cpos * (pi/180)),3),
                         " (S", 2:nc, ii, ");"), 
                  "")
    }
    
    inpM <- c(inpM,
              
    "\t\t\t!FACTOR (CO-)VARIANCES:", 
    
    "",
    
    paste0("\t\t\t", vc[1], " - ", vc[nc], "@0;"),
    
    paste0("\t\t\tL_",vc[1], " - L_", vc[nc], "@0;"),
    
    paste0("\t\t\tU_",vc[1], " - U_", vc[nc], "@1;"),
    
    paste0("\t\t\tLEV - SIN", fcomp, "@1;"),
    
    paste0("\t\t\tU_",vc[1], " - U_", vc[nc], " WITH ",
           "U_", vc[1], " - U_", vc[nc],"@0;"),
    
    paste0("\t\t\tU_",vc[1], " - U_", vc[nc], " WITH ",
           "LEV - SIN", fcomp, "@0;"), 
    
    paste0("\t\t\tLEV - SIN", fcomp, " WITH LEV - SIN", fcomp, "@0;"),
    
    "")
    
    if(ne > 0){
        inpM <- c(inpM, 
            "\t\t\t!MEASUREMENT PART FOR COVARIATES", "")
        
        for(ii in 1:ne){
            
            inpM <- c(inpM,
            
            paste0("\t\t\tL_", ve[ii], " BY ", ve[ii], "*0.8;" ),
            
            paste0("\t\t\t", ve[ii], "@0;" ),
            
            paste0("\t\t\tL_", ve[ii], "@1;" ),
            
            "")
        }
        
        if(ne > 1){
            for(ii in 1:nrow(svalCovariateCor)){
                inpM <- c(inpM,
                          paste0("\t\t\tL_EV", svalCovariateCor[ii,"v1"],
                                 " WITH L_EV", svalCovariateCor[ii,"v2"],
                                 "*", svalCovariateCor[ii,"sval"], ";")
                          )
            }
        }
        
        inpM <- c(inpM, "", 
                  
            "\t\t\t!RELATIONSHIP(S) BETWEEN COMMON AND UNIQUE FACTORS WITH COVARIATE(S)", 
            
            "\t\t\t!COMMON FACTORS:", ""
            
            )
        
        for (jj in 1:ne) {
            
            inpM <- c(inpM, paste0("\t\t\tLEV  WITH L_", ve[jj], "*.3 (rl0", letters[jj],");"))
                      
            for(ii in 1:fcomp){
            inpM <- c(inpM,
                
                paste0("\t\t\tCOS", ii, " WITH L_", ve[jj], "*.3 (rc", ii, letters[jj], ");"),
                
                paste0("\t\t\tSIN", ii, " WITH L_", ve[jj], "*.3 (rs", ii, letters[jj], ");")
            )
            }
            
            inpM <- c(inpM, "")
        }
        
        inpM <- c(inpM, "\t\t\t!UNIQUE FACTORS", "")
        
        svalTemp <- rep(c(.05, -.05), times=nc/2)
        
        for (ii in 1:ne) {
            inpM <- c(
                inpM,
                paste0("\t\t\tU_", vc, " WITH L_", ve[ii], "*", svalTemp, "(ru", 1:nc, letters[ii], ");"),
                "")
        }
    }
    
    inpM <- c(inpM,
    
    "\t\t\tMODEL CONSTRAINT:",
    
    "\t\t\t!NEW PARAMETERS PERTAINING TO THE CORRELATION FUNCTION:",
    
    "",
    
    paste0("\t\t\tNEW(t", 2:nc, "*", round(start.cpos * (pi/180),3), ");"),
    
    "",
    
    paste0("\t\t\tt", 2:nc, " < ", 2*pi, ";"),
    
    "",
    
    paste0("\t\t\tt", 2:nc, " > 0;"),
    
    "",
   
    "\t\t\t!CONSTRAINING BETA-PARAMETERS",
    
    "",
    
    paste0("\t\t\tb0 = sqrt(1 - (", paste0("b", 1:fcomp, "**2", collapse=" + "), "));"),
    
    "",
    
    paste0("\t\t\tb", 1:fcomp, " > 0;"),
    
    "")
    
    for(ii in 1:fcomp){
        inpM <- c(inpM,
                  paste0("\t\t\tc", 2:nc, ii, " = b", ii, "*cos(",ii, "*t", 2:nc, "); \t",
                         "s", 2:nc, ii, " = b", ii, "*sin(",ii, "*t", 2:nc, ");"), "")
    }
    
    # CONSTRAINTS FOR RELATIONSHIPS WITH COVARIATES
    if(ne > 0){
        
        for (ii in 1:ne) {
            inpM <- c(
                inpM,
                paste0("\t\t\tNEW(g0", letters[ii], "*.5);"),
                paste0("\t\t\tNEW(g1", letters[ii], "*.5);"),
                paste0("\t\t\tNEW(d", letters[ii], "*", start.epos[ii], ");"), 
                "", 
                paste0("\t\t\tg1", letters[ii], " > 0;"),
                paste0("\t\t\td", letters[ii], " < ", 2*pi, ";"),
                paste0("\t\t\td", letters[ii], " > 0;"),
                ""
            )
        }
        
        for (ii in 1:ne) {
            
            inpM <- c(inpM, paste0("\t\t\trl0", letters[ii], " = g0", letters[ii], "*b0;"))
            
            for (jj in 1:fcomp) {
                inpM <- c(
                    inpM,
                    paste0("\t\t\trc", jj, letters[ii], " = g1", letters[ii], "*b", jj,"*cos(", jj, "*d", letters[ii],");"),
                    paste0("\t\t\trs", jj, letters[ii], " = g1", letters[ii], "*b", jj,"*sin(", jj, "*d", letters[ii],");")
                )
            }
            
            inpM <- c(inpM, "")
        }
        
        # CONSTRAINTS FOR BETA_0
        for(ii in 1:ne){
            inpM <- c(inpM, "\t\t\t0 = ")
            for(jj in 1:nc){
                if(jj < nc){
                    inpM <- c(inpM, paste0("\t\t\t(b0*ru", jj, letters[ii], ")/U", jj, "+"))
                } else if(jj == nc) {
                    inpM <- c(inpM, paste0("\t\t\t(b0*ru", jj, letters[ii], ")/U", jj, ";"), "")
                }
            }
        }
        
        # CONSTRAINTS FOR OTHER BETA(S)
        for(ii in 1:ne){
            
            inpM <- c(inpM, "\t\t\t0 = ")
            
            for(jj in 1:nc){
                
                    if(jj==1){
                        
                        inpM <- c(inpM, paste0(c("\t\t\t((", 
                                          paste0("b", 1:fcomp, "**2*cos(", 1:fcomp, "*(0 -d", letters[ii],"))",collapse="+\n\t\t\t" ),
                                          paste0(")*ru", jj, letters[ii], ")/U", jj, "+")), collapse=""))
                        
                    } else if(jj > 1 & jj < nc) {
                        
                        inpM <- c(inpM, paste0(c("\t\t\t((", 
                                                 paste0("b", 1:fcomp, "**2*cos(", 1:fcomp, "*(t", jj, "-d", letters[ii],"))",collapse="+\n\t\t\t" ),
                                                 paste0(")*ru", jj, letters[ii], ")/U", jj, "+")), collapse=""))
                        
                    } else if(jj == nc){
                        
                        inpM <- c(inpM, paste0(c("\t\t\t((", 
                                                 paste0("b", 1:fcomp, "**2*cos(", 1:fcomp, "*(t", jj, "-d", letters[ii],"))",collapse="+\n\t\t\t" ),
                                                 paste0(")*ru", jj, letters[ii], ")/U", jj, ";")), collapse=""), "")
                        
                }
            }
        }
        
        for(ii in 1:ne){
            
            inpM <- c(inpM, "\t\t\t0 = ")
            
            for(jj in 1:nc){
                
                if(jj==1){
                    
                    inpM <- c(inpM, paste0(c("\t\t\t((", 
                                             paste0(1:fcomp, "*b", 1:fcomp, "**2*sin(", 1:fcomp, "*(d", letters[ii],"-0 ))",collapse="+\n\t\t\t" ),
                                             paste0(")*ru", jj, letters[ii], ")/U", jj, "+")), collapse=""))
                    
                } else if(jj > 1 & jj < nc) {
                    
                    inpM <- c(inpM, paste0(c("\t\t\t((", 
                                             paste0(1:fcomp, "*b", 1:fcomp, "**2*sin(", 1:fcomp, "*(d", letters[ii], "-t", jj,"))", collapse="+\n\t\t\t" ),
                                             paste0(")*ru", jj, letters[ii], ")/U", jj, "+")), collapse=""))
                    
                } else if(jj == nc){
                    
                    inpM <- c(inpM, paste0(c("\t\t\t((", 
                                             paste0(1:fcomp, "*b", 1:fcomp, "**2*sin(", 1:fcomp, "*(d", letters[ii], "-t", jj,"))", collapse="+\n\t\t\t" ),
                                             paste0(")*ru", jj, letters[ii], ")/U", jj, ";")), collapse=""), "")
                    
                }
            }
        }
    }
    
    # OPTIONAL CONSTRAINTS
    if (eq.scaling == TRUE) {
        inpM <- c(
            inpM,
            "\t\t\t!EQUAL SCALING PARAMETERS",
            paste0("\t\t\tL1 = L", 2:nc, ";"),
            ""
        )
    }
    
    if (eq.unique == TRUE) {
        inpM <- c(inpM,
                  "\t\t\t!EQUAL UNIQUE VARIANCES",
                  paste0("\t\t\tU1 = U", 2:nc, ";"),
                  "")
    }
    
    if (pos.fixed == TRUE) {
        inpM <- c(inpM,
                  "\t\t\t!FIXING INDICATOR SCALE POSITIONS",
                  paste0("\t\t\tt", 2:nc, " = ", round(start.cpos * (pi / 180), 3), ";"),
                  "")
    }
    
    # CONSTRAINTS FOR INTERVAL ESTIMATES
    if(ne > 0){
        
        inpM <- c(inpM, "\t\t\t!PHI PARAMETERS FOR CONFIDENCE ELLIPSES AND WALD TESTS", "")
        
        for(ii in 1:ne){
            inpM <- c(inpM,
                      paste0("\t\t\tNEW(phi_c", letters[ii], "*.1 phi_s", letters[ii], "*.1);"),
                      paste0("\t\t\tphi_c", letters[ii], " = g1", letters[ii], "*cos(d", letters[ii],");"),
                      paste0("\t\t\tphi_s", letters[ii], " = g1", letters[ii], "*sin(d", letters[ii],");"),
                      ""
                      )
        }
    }
    
    # WALD TESTS 
    wald.test <- match.arg(wald.test)

    if(!wald.test=="none"){
        
        if(length(wald.vars) < 2 | !all(wald.vars %in% evars)){
          
          stop("Argument 'wald.vars' must be a subset of argument 'evars' with at least 2 valid entries")
          
        }

        inpM <- c(inpM, 
                  paste0("\t\t\t!WALD TEST FOR CALL: ", toupper(wald.test)),
                  "", "\t\t\tMODEL TEST:", "")
    
        idx <- sort(pmatch(wald.vars, evars))
    
        nwv <- length(idx)
    
        if(wald.test == "coordinate"){
            for(ii in 2:nwv){
                inpM <- c(inpM,
                          paste0("\t\t\tphi_c", letters[idx[1]], " = phi_c", letters[idx[ii]], ";"),
                          paste0("\t\t\tphi_s", letters[idx[1]], " = phi_s", letters[idx[ii]], ";")
                        )
            }
        } else if(wald.test == "sensitivity"){
            for(ii in 2:nwv){
              inpM <- c(inpM,
                        paste0("\t\t\tg1", letters[idx[1]], " = g1", letters[idx[ii]], ";")
              )
            }
        } else if(wald.test == "zero_sens"){
            for(ii in 1:nwv){
              inpM <- c(inpM,
                        paste0("\t\t\t0 = g1", letters[ii], ";")
              )
            }
        } else if(wald.test == "displacement"){
          for(ii in 2:(nwv)){
            inpM <- c(inpM,
                      paste0("\t\t\t0 = 1 - cos(d", letters[idx[1]], " - d", letters[idx[ii]], ");")
            )
          }
        }
    }

    inpM <- c(inpM, "", 
              "Output:\t\tTECH1 TECH3;")
    
    # OUTPUT
    writeLines(inpM, con=paste0(file,".inp"))
    
    message(paste0("Syntax created as ", file,".inp"))
    
    # RUN MODEL AND CHECK
    if(run==TRUE){
        runModels(paste0(file, ".inp"), showOutput=TRUE)
      
        checkout <- readLines(paste0(file, ".out"))
        
        if(any(grepl("NO CONVERGENCE", checkout))){
          message("\nThe Mplus algorithm did not converge.\nTry adjusting the starting values pertaining to the covariates' positions.")
        } else if (any(grepl("NON-POSITIVE", checkout))){
          message("\nMplus probably did not find an appropriate solution.\nTry adjusting the starting values pertaining to the covariates' positions.")
        } else {
          message("\nThe Mplus algorithm has converged successfully.\nThe solution appears to be appropriate.")
        }
    }
}



