do.kin.blup <- function(pheno_data, geno_data, trait, geno_dist, env_to_mask, prop, cycles, multi_env = NULL){ 
    
    perct <- prop/100;
    work_pheno <- pheno_data
    work_geno <- geno_data
    num <- 100 - prop
      
    if(prop > 0 && prop < 100){                
        cat(blue("Adding", num, "percent of the val pop to TP\n"))           
        traits <- 3; cycles <- cycles              
        pred_ability_df <- matrix(nrow = cycles, ncol = traits)        
        for(r in 1:cycles){       
            pheno_NA <- work_pheno
            val_pop <- work_pheno %>% filter(Env == env_to_mask)
            val_NA <- sample(rownames(pheno_NA[pheno_NA$Env == env_to_mask, ]), round(nrow(val_pop)*perct))
            pheno_NA[val_NA, trait] <- NA
            fit <- kin.blup(pheno_NA,K=geno_dist,GAUSS=TRUE,geno="ID",pheno=trait, fixed = multi_env)
            temp1 <- as.data.frame(fit$pred); temp2 <- cbind(rownames(temp1), temp1)            
            predicted_trait <- paste0("predicted_", trait); colnames(temp2) <- c("ID", predicted_trait)             
            temp3 <- suppressMessages(temp2 %>% inner_join(val_pop)); pheno_NA_v <- pheno_NA %>% filter(Env == env_to_mask)            
            pheno_NA_v <- pheno_NA_v[,c("ID", trait)]; colnames(pheno_NA_v) <- c("ID", "Masked_trait")            
            temp4 <- suppressMessages(temp3 %>% inner_join(pheno_NA_v)) 
            temp5 <- temp4 %>% filter(is.na(Masked_trait))
            pred_ability_df[r, 1] <- cor(temp4[, predicted_trait], temp4[, trait], use = "complete.obs")
            pred_ability_df <- as.data.frame(pred_ability_df); colnames(pred_ability_df) <- paste0("Pearson.", trait)                    
            out <- pred_ability_df                   
        }
        return(out)
    }else if(prop == 100){        
        cat(blue("Adding", num, "percent of the val pop to TP\n"))          
        pheno_NA <- work_pheno
        val_pop <- work_pheno %>% filter(Env == env_to_mask)
        val_NA <- sample(rownames(pheno_NA[pheno_NA$Env == env_to_mask, ]), round(nrow(val_pop)*perct))
        pheno_NA[val_NA, trait] <- NA
        fit <- kin.blup(pheno_NA,K=geno_dist,GAUSS=TRUE,geno="ID",pheno=trait, fixed = multi_env)
        temp1 <- as.data.frame(fit$pred); temp2 <- cbind(rownames(temp1), temp1)        
        predicted_trait <- paste0("predicted_", trait); colnames(temp2) <- c("ID", predicted_trait)         
        temp3 <- suppressMessages(temp2 %>% inner_join(val_pop)); pheno_NA_v <- pheno_NA %>% filter(Env == env_to_mask)        
        pheno_NA_v <- pheno_NA_v[,c("ID", trait)]; colnames(pheno_NA_v) <- c("ID", "Masked_trait")        
        temp4 <- suppressMessages(temp3 %>% inner_join(pheno_NA_v)) 
        pred_ability_df <- cor(temp4[, predicted_trait], temp4[, trait], use = "complete.obs")
        pred_ability_df <- as.data.frame(pred_ability_df); colnames(pred_ability_df) <- paste0("Pearson.", trait)                
        out <- pred_ability_df        
        return(out)
    }else if(prop == 0){        
        cat(blue("Adding", num, "percent of the val pop to TP\n"))         
        pheno_NA <- work_pheno
        val_pop <- work_pheno %>% filter(Env == env_to_mask)
        val_NA <- sample(rownames(pheno_NA[pheno_NA$Env == env_to_mask, ]), round(nrow(val_pop)*perct))
        pheno_NA[val_NA, trait] <- NA
        fit <- kin.blup(pheno_NA,K=geno_dist,GAUSS=TRUE,geno="ID",pheno=trait, fixed = multi_env)
        temp1 <- as.data.frame(fit$pred); temp2 <- cbind(rownames(temp1), temp1)        
        predicted_trait <- paste0("predicted_", trait); colnames(temp2) <- c("ID", predicted_trait)         
        temp3 <- suppressMessages(temp2 %>% inner_join(val_pop)); pheno_NA_v <- pheno_NA %>% filter(Env == env_to_mask)        
        pheno_NA_v <- pheno_NA_v[,c("ID", trait)]; colnames(pheno_NA_v) <- c("ID", "Masked_trait")        
        temp4 <- suppressMessages(temp3 %>% inner_join(pheno_NA_v)) 
        pred_ability_df <- cor(temp4[, predicted_trait], temp4[, trait], use = "complete.obs")
        pred_ability_df <- as.data.frame(pred_ability_df); colnames(pred_ability_df) <- paste0("Pearson.", trait)                
        out <- pred_ability_df        
        return(out)        
    }
  }      
