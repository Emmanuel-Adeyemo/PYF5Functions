do.mixed.solve <- function(pheno_for_train, geno_for_train, fixed_effect_matrix, pheno_for_val, geno_for_val, trait){
    
    # Deal with input
    pheno_train_trait <- pheno_for_train[, trait]
    pheno_train_trait <- as.matrix(pheno_train_trait)
    work_geno_train <- as.matrix(geno_for_train[,2:ncol(geno_for_train)])
    work_geno_val <- geno_for_val
    rownames(work_geno_val) <- work_geno_val[,1]
    work_geno_val <- work_geno_val[,-1]  
    work_geno_val <- as.matrix(work_geno_val)
    
    # Solve the mixed model
    solve_out <- mixed.solve(y = pheno_train_trait, Z = work_geno_train, X = fixed_effect_matrix, method = "REML")
    marker_effects <- solve_out$u

    # Predict genotypic values
    geno_values <- work_geno_val %*% marker_effects
    geno_values <- geno_values[,1] + solve_out$beta[1]
    geno_values <- as.data.frame(geno_values)
    geno_values <- cbind(ID = rownames(geno_values), geno_values)
    predicted_val <- suppressMessages(geno_values %>% inner_join(pheno_for_val))

    # find correlations between the predicted and observed values in the validation population
    pred_accuracy <- cor(predicted_val$geno_values, predicted_val[,trait], use = "complete.obs")
    
    out_pred <- predicted_val %>% dplyr::select(ID, geno_values)
    genoo_value_trait <- paste0("geno_values_", trait)
    colnames(out_pred) <- c("ID", genoo_value_trait)
    out <- list(pred_accuracy, out_pred)
    # Return the data
    return(out)    
} 
