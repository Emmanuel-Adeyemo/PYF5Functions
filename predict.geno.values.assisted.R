predict.geno.values.assisted <- function(full_geno, pheno_train, geno_train, pheno_val, geno_val, prop, cycles) { # function hardcorded for Traits

    # libraries to load
    require(dplyr)
    require(rrBLUP)
    require(purrr)
    require(gsubfn)
    
    # checks if this is the appropiate fuction to use
    if(prop < 1 || prop > 99){
        stop(cat(blue("This function cannot add", prop, "percent of the val pop to TP, use predict.geno.values() instead.\n")))
    }
    start <- Sys.time()
    # Deal with input
    perct <- prop/100;
    num <- round(nrow(geno_val)*perct)
    traits <- 3      # Hardcoded
    
    pred_ability_df <- matrix(nrow = cycles, ncol = traits) 
    twt_list <- list(); vsk_list <- list(); dis_list <- list()   
    colnames(pred_ability_df) <- c("Pearson.twt", "Pearson.vsk", "Pearson.dis" )
    
    # Loop starts here
    for(r in 1:cycles){
        
    samp_train <- sample(nrow(pheno_train), num)
    new_train <- train_17_pheno[-samp_train,]
    new_train_ge <- suppressMessages(new_train %>% dplyr::select(ID) %>% inner_join(full_geno))
    
    samp_val <- sample(nrow(geno_val), num)
    new_val_geno <- geno_val[-samp_val,]
    take_out_val_geno <- geno_val[samp_val,]
    val_id_to_add <- take_out_val_geno %>% dplyr::select(ID)
    val_id_to_add_pheno <- suppressMessages(val_id_to_add %>% inner_join(pheno_val))
    
    new_train_pheno <- rbind(new_train, val_id_to_add_pheno)
    new_train_geno <- rbind(new_train_ge, take_out_val_geno)
       
    XF <- as.matrix(model.matrix(ID~Env, data = new_train_pheno))
    
    list[acc.twt, pgv.twt] <- do.mixed.solve(new_train_pheno, new_train_geno, XF, pheno_val, new_val_geno, "Twt")  # Trait 1
    list[acc.vsk, pgv.vsk] <- do.mixed.solve(new_train_pheno, new_train_geno, XF, pheno_val, new_val_geno, "VSK")  # Trait 2
    list[acc.dis, pgv.dis] <- do.mixed.solve(new_train_pheno, new_train_geno, XF, pheno_val, new_val_geno, "DIS")  # Trait 3
    
    pred_ability_df[r, 1] <-  acc.twt; pred_ability_df[r, 2] <-  acc.vsk; pred_ability_df[r, 3] <-  acc.dis 
    twt_list[[length(twt_list) + 1]] <- pgv.twt; vsk_list[[length(vsk_list) + 1]] <- pgv.vsk; dis_list[[length(dis_list) + 1]] <- pgv.dis
      
    } # Loop ends here
    twt_reduce <- twt_list %>% reduce(full_join, by = "ID")
    vsk_reduce <- vsk_list %>% reduce(full_join, by = "ID")
    dis_reduce <- dis_list %>% reduce(full_join, by = "ID")
    
    out <- list(pred_ability_df, twt_reduce, vsk_reduce, dis_reduce)
    end <- Sys.time()
    print(end - start)
    return(out)
} 
