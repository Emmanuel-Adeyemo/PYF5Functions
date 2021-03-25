predict.validation.random.multi.trait <- function(full_geno, pheno_train, geno_train, pheno_val, tpSize, cycles) { # function hardcorded for Traits

    require(dplyr)
    require(MTM)
    require(reshape2)
    require(MCMCpack)
    
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
        
        vpSet <- sample_n(pheno_val, vpSize)
        vpSetGeno <- vpSet %>% dplyr::select(ID) %>% inner_join(full_geno,  by = "ID")
        remaining_pheno <- pheno_val %>% anti_join(vpSet, by = "ID")
        
        tpSet <- sample_n(remaining_pheno, tpSize)
        tpSetGeno <- tpSet %>% dplyr::select(ID) %>% inner_join(full_geno,  by = "ID")
        tpLargeSet <- rbind(pheno_train, tpSet)
        tpLargeSetGeno <- tpLargeSet %>% dplyr::select(ID) %>% inner_join(full_geno,  by = "ID")
                
        smallSetGeno <- rbind(tpSetGeno, vpSetGeno)
        smallSetGeno <- smallSetGeno %>% arrange(ID)    
        largeSetGeno <- rbind(tpLargeSetGeno, vpSetGeno)
        largeSetGeno <- largeSetGeno %>% arrange(ID)

        vpSetNA <- vpSet
        vpSetNA$Twt <- vpSetNA$VSK <- vpSetNA$DIS <- NA

        smallPheno <- rbind(tpSet, vpSet)
        smallPheno <- smallPheno %>% arrange(ID)

        smallPhenoNA <- rbind(tpSet, vpSetNA)
        smallPhenoNA <- smallPhenoNA %>% arrange(ID)

        largePheno <- rbind(tpLargeSet, vpSet)
        largePheno <- largePheno %>% arrange(ID)

        largePhenoNA <- rbind(tpLargeSet, vpSetNA)
        largePhenoNA <- largePhenoNA %>% arrange(ID)
        
        accuracyLarge[r, ] <- multi.trait.prediction(largeSetGeno, largePheno, largePhenoNA)
        accuracySmall[r, ] <- multi.trait.prediction(smallSetGeno, smallPheno, smallPhenoNA)
      
    } 
        
    out <- list(accuracySmall, accuracyLarge)
    end <- Sys.time()
    print(end - start)
    return(out)
} 
