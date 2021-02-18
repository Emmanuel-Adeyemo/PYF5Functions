select.by.genetic.algorithm <- function(validation_genotpe, num_to_select, method){
    
    #'@ method: options are "AOPT", "DOPT", "PEVMEAN", "CDMEAN" and so on. See STPGA manual for details
    #'@ output is a dataframe of optimized lines ie TP
    
    rownames(validation_genotpe) <- validation_genotpe[, 1]
    validation_genotpe <- validation_genotpe[,-1]
    
    M <- validation_genotpe
    Kmat<-cov(t(M))
    Kmat<-Kmat/mean(diag(Kmat))
    dim(Kmat)
    svdG<-svd(Kmat, nu=40,nv=40)
    PCAs<-Kmat%*%svdG$v
    rownames(PCAs)<-rownames(Kmat)

    ListTrain1<-GenAlgForSubsetSelectionNoTest(P=PCAs,ntoselect=num_to_select, InitPop=NULL,
                                                npop=500, nelite=10, mutprob=.5, mutintensity = 1,
                                                niterations=800,minitbefstop=20, tabu=F,tabumemsize = 0,
                                                plotiters=F, lambda=1e-9,errorstat=method, mc.cores=1)
    X <- ListTrain1[[1]]
    out <- as.data.frame(X)
    colnames(out) <- "ID"
    
    return(out)
    
}
