multi.trait.prediction <- function(genoSet, phenoSet, phenoSetNA){
    
    rownames(genoSet) <- genoSet[,1]
    genoSet <- genoSet[,-1]
    genoSet <- genoSet * 2 -1 
    
    nTraits <- 3
    
    pred_ability = matrix(nrow=1, ncol=nTraits)
    
    geno <- as.matrix(genoSet)
        
    phenoSplit <- phenoSetNA[,c(1,2:4)]
    phenoSplit <- as.matrix(phenoSplit[,2:ncol(phenoSplit)])
    
    phenoSplitDiag <- diag(1, nrow(phenoSplit), nrow(phenoSplit))
    
    G <- tcrossprod(geno)/ncol(geno)
    
    multimodel <- MTM(Y = phenoSplit,
                K = list(list(K = G, COV = list(type = 'UN', df0 = nTraits, S0 = diag(nTraits))), 
                         list(K = phenoSplitDiag, COV = list(type = 'UN', df0 = nTraits, S0 = diag(nTraits)))),
                resCov = list(type = 'DIAG', S0 = rep(1, nTraits), df0 = rep(1, nTraits)), 
                nIter = 100, burnIn = 25, thin = 3, saveAt = 'ex1')
    
    tmp1 <- as.data.frame(multimodel$YHat)
    tmp2 <- as.data.frame(cbind(phenoSet[,2], phenoSplit[,1], tmp1[,1]))
    tmp3 <- tmp2[rowSums(is.na(tmp2)) > 0,]
    pred_ability[1,1] <- cor(tmp3$V1, tmp3$V3, use="complete.obs")

    tmp2 <- as.data.frame(cbind(phenoSet[,3], phenoSplit[,2], tmp1[,2]))
    tmp3 <- tmp2[rowSums(is.na(tmp2)) > 0,]
    pred_ability[1,2] <- cor(tmp3$V1, tmp3$V3, use="complete.obs")

    tmp2 <- as.data.frame(cbind(phenoSet[,4], phenoSplit[,3], tmp1[,3]))
    tmp3 <- tmp2[rowSums(is.na(tmp2)) > 0,]
    pred_ability[1,3] <- cor(tmp3$V1, tmp3$V3, use="complete.obs")
    
    return(pred_ability)
    
}
