predict.geno.values <- function(pheno_train, geno_train, pheno_val, geno_val, trait) {

    # load required libraries
    require(dplyr)
    require(rrBLUP)

    # Deal with input
    pheno_train_trait <- pheno_train[, trait]
    pheno_train_trait <- as.matrix(pheno_train_trait)
    work_geno_train <- as.matrix(geno_train[,2:ncol(geno_train)])
    work_geno_val <- geno_val
    rownames(work_geno_val) <- work_geno_val[,1]
    work_geno_val <- work_geno_val[,-1]  
    work_geno_val <- as.matrix(work_geno_val)
    
    XF <- as.matrix(model.matrix(ID~Env, data = pheno_train))

    # Solve the mixed model
    solve_out <- mixed.solve(y = pheno_train_trait, Z = work_geno_train, X = XF, method = "REML")
    marker_effects <- solve_out$u

    # Predict genotypic values
    geno_values <- work_geno_val %*% marker_effects

    row.names(geno_values) <- row.names(geno_val)

    # find correlations between the predicted and observed values in the validation population
    pred_accuracy <- cor(geno_values, pheno_val[,trait], use = "complete.obs")


    # Return the data
    return(pred_accuracy)

} 
