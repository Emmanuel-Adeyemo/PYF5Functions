predict.validation.family.assisted <- function(full_geno, pheno_train, geno_train, pheno_val, val_pedigree, geno_val, tpSize, cycles) { # function hardcorded for Traits

    require(dplyr)
    require(rrBLUP)
    require(purrr)
    require(gsubfn)
    
    if(tpSize > nrow(pheno_val)){
        stop(cat(blue("TP size of", tpSize, "is too large. Size must be lower than VP size  \n")))
    }
    start <- Sys.time()
    # Deal with input
    
    maxtpSize <- 200
    vpSize <- nrow(pheno_val) - maxtpSize
    
    traits <- 3      # Hardcoded for the three traits in this study
    
    accuracyLarge <- matrix(nrow = cycles, ncol = traits) 
    accuracySmall <- matrix(nrow = cycles, ncol = traits)
    
    colnames(accuracyLarge) <- colnames(accuracySmall) <- c("Pearson.twt", "Pearson.vsk", "Pearson.dis" )
    
    # Loop starts here
    for(r in 1:cycles){
        
        tpSet <- sample.validation.by.pedigree(pheno_val, val_pedigree, tpSize)
        tpSet <- tpSet %>% dplyr::select(-Pedigree)
        tpSetGeno <- tpSet %>% dplyr::select(ID) %>% inner_join(full_geno,  by = "ID") # get geno for selected tp
        
        tpLargeSet <- rbind(pheno_train, tpSet) # larger TP
        tpLargeSetGeno <- tpLargeSet %>% dplyr::select(ID) %>% inner_join(full_geno,  by = "ID") # get geno for larger TP
        
        remaining_pheno <- pheno_val %>% anti_join(tpSet, by = "ID") # remaining candidates
        vpSet <- sample_n(remaining_pheno, vpSize) # randomly select a fixed size VP
        vpSetGeno <- vpSet %>% dplyr::select(ID) %>% inner_join(full_geno,  by = "ID") # get geno for the random VP 

        accuracySmall[r, 1] <- do.mixed.solve(tpSet, tpSetGeno, vpSet, vpSetGeno, "Twt")
        accuracySmall[r, 2] <- do.mixed.solve(tpSet, tpSetGeno, vpSet, vpSetGeno, "VSK")
        accuracySmall[r, 3] <- do.mixed.solve(tpSet, tpSetGeno, vpSet, vpSetGeno, "DIS")
        
        accuracyLarge[r, 1] <- do.mixed.solve(tpLargeSet, tpLargeSetGeno, vpSet, vpSetGeno, "Twt")
        accuracyLarge[r, 2] <- do.mixed.solve(tpLargeSet, tpLargeSetGeno, vpSet, vpSetGeno, "VSK")
        accuracyLarge[r, 3] <- do.mixed.solve(tpLargeSet, tpLargeSetGeno, vpSet, vpSetGeno, "DIS")   
      
    } 
    out <- list(accuracySmall, accuracyLarge)
    end <- Sys.time()
    print(end - start)
    return(out)
} 
