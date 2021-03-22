predict.validation.genetic.assisted <- function(full_geno, pheno_train, geno_train, pheno_val, geno_val, tpSize, cycles, method) { # function hardcorded for Traits

    require(dplyr)
    require(rrBLUP)
    require(purrr)
    require(gsubfn)
    require(STPGA)
    
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
    
    twt_small <- vsk_small <- dis_small <- list()   
    twt_large <- vsk_large <- dis_large <- list() 
    
    # Loop starts here
    for(r in 1:cycles){
        
        tpSetList <- select.by.genetic.algorithm(geno_val, tpSize, method)  # select TP from candidates
        
        tpSet <- tpSetList %>% inner_join(pheno_val, by = "ID") # get pheno for small TP
        tpSetGeno <- tpSet %>% dplyr::select(ID) %>% inner_join(full_geno,  by = "ID") # get geno for selected tp
        
        tpLargeSet <- rbind(pheno_train, tpSet) # larger TP
        tpLargeSetGeno <- tpLargeSet %>% dplyr::select(ID) %>% inner_join(full_geno,  by = "ID") # get geno for larger TP
        
        remaining_pheno <- pheno_val %>% anti_join(tpSet, by = "ID") # remaining candidates
        
        vpSet <- sample_n(remaining_pheno, vpSize) # randomly select a fixed size VP
        vpSetGeno <- vpSet %>% dplyr::select(ID) %>% inner_join(full_geno,  by = "ID") # get geno for the random VP      
        
        # small TP
        list[accSmall.twt, pgvSmall.twt] <- do.mixed.solve(tpSet, tpSetGeno, vpSet, vpSetGeno, "Twt") # Trait 1
        list[accSmall.vsk, pgvSmall.vsk] <- do.mixed.solve(tpSet, tpSetGeno, vpSet, vpSetGeno, "VSK") # Trait 2
        list[accSmall.dis, pgvSmall.dis] <- do.mixed.solve(tpSet, tpSetGeno, vpSet, vpSetGeno, "DIS") # Trait 3

        accuracySmall[r, 1] <-  accSmall.twt; accuracySmall[r, 2] <-  accSmall.vsk; accuracySmall[r, 3] <-  accSmall.dis 
        
        twt_small[[length(twt_small) + 1]] <- pgvSmall.twt; vsk_small[[length(vsk_small) + 1]] <- pgvSmall.vsk; dis_small[[length(dis_small) + 1]] <- pgvSmall.dis
        
        # large TP
        list[accLarge.twt, pgvLarge.twt] <- do.mixed.solve(tpLargeSet, tpLargeSetGeno, vpSet, vpSetGeno, "Twt") # Trait 1
        list[accLarge.vsk, pgvLarge.vsk] <- do.mixed.solve(tpLargeSet, tpLargeSetGeno, vpSet, vpSetGeno, "VSK") # Trait 2
        list[accLarge.dis, pgvLarge.dis] <- do.mixed.solve(tpLargeSet, tpLargeSetGeno, vpSet, vpSetGeno, "DIS") # Trait 3

        accuracyLarge[r, 1] <-  accLarge.twt; accuracyLarge[r, 2] <-  accLarge.vsk; accuracyLarge[r, 3] <-  accLarge.dis 
        
        twt_large[[length(twt_large) + 1]] <- pgvLarge.twt; vsk_large[[length(vsk_large) + 1]] <- pgvLarge.vsk; dis_large[[length(dis_large) + 1]] <- pgvLarge.dis
      
    } 
    twtRedSmall <- twt_small %>% reduce(full_join, by = "ID"); vskRedSmall <- vsk_small %>% reduce(full_join, by = "ID"); disRedSmall <- dis_small %>% reduce(full_join, by = "ID")
    twtRedLarge <- twt_large %>% reduce(full_join, by = "ID"); vskRedLarge <- vsk_large %>% reduce(full_join, by = "ID"); disRedLarge <- dis_large %>% reduce(full_join, by = "ID")
    
    RedSmall <- list(twtRedSmall, vskRedSmall, disRedSmall)
    RedLarge <- list(twtRedLarge, vskRedLarge, disRedLarge)
    
    out <- list(accuracySmall, accuracyLarge, RedSmall, RedLarge)
    end <- Sys.time()
    print(end - start)
    return(out)
} 
