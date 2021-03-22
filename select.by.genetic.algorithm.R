select.by.genetic.algorithm <- function(validation_genotpe, num_to_select, method){
    
    #'@ method: options are "AOPT", "DOPT", "PEVMEAN", "CDMEAN", "CDMEANMM (for multitrait)"
    
    if(method == "CDMEANMM"){        
        Lambda <- 1
    }else{
        Lambda <- 1e-9
    }
    
    rownames(validation_genotpe) <- validation_genotpe[, 1]
    validation_genotpe <- validation_genotpe[,-1]
    
    M <- validation_genotpe
    Kmat<-cov(t(M))
    Kmat<-Kmat/mean(diag(Kmat))
    dim(Kmat)
    svdG<-svd(Kmat, nu=40,nv=40)
    PCAs<-Kmat%*%svdG$v
    rownames(PCAs)<-rownames(Kmat)
    
    if(method == "CDMEAN" | method == "CDMEANMM"){
        Psolve <- solve(Kmat+1e-6*diag(ncol(Kmat)))
    }else{
        Psolve <- PCAs
    }
    ListTrain1<-GenAlgForSubsetSelectionNoTest(P=Psolve,ntoselect=num_to_select, InitPop=NULL,
                                                npop=200, nelite=10, mutprob=.5, mutintensity = 1,
                                                niterations=500,minitbefstop=20, tabu=F,tabumemsize = 0,
                                                plotiters=F, lambda=Lambda,errorstat=method, mc.cores=1) 
    X <- ListTrain1[[1]]
    out <- as.data.frame(X)
    colnames(out) <- "ID"
    
    return(out)
    
}
