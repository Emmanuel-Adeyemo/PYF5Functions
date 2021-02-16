sample.validation.by.pedigree <- function(val_pheno, val_pedigree_data, prop_num){
    
    val_data_fr <- val_pedigree_data %>% inner_join(val_pheno, by = "ID")
    in_order <- val_data_fr %>% count(Pedigree) %>% arrange(desc(n))
    maxn <- max(in_order$n)
    selected_df_by_pedigree <- as.data.frame(matrix(ncol = 7))
    colnames(selected_df_by_pedigree) <- c("ID", "Pedigree", "Twt", "VSK", "DIS", "Env", "GID")
    selected_df_by_pedigree <- selected_df_by_pedigree[-1,]
    selected_df_by_pedigree$ID <- as.factor(selected_df_by_pedigree$ID)
    nround <- 1
    while(nround < maxn){

        rows_remaining <- prop_num - nrow(selected_df_by_pedigree)
        peds_more_than_round <- in_order %>% filter(n > nround)
        
        if(rows_remaining == 0){
            break
        }
        
        if(rows_remaining >= nrow(peds_more_than_round)){
            selected_ped <- in_order[1:nrow(peds_more_than_round), 1]
        }else if(rows_remaining < nrow(peds_more_than_round)){
            selected_ped <- in_order[1:rows_remaining, 1]
        }    
            
        for(i in selected_ped){
            val_data_remaining <- val_data_fr %>% anti_join(selected_df_by_pedigree, by = "ID")
            selected_rows <- sample_n(subset(val_data_remaining, Pedigree == i), 1)
            selected_df_by_pedigree <- rbind(selected_df_by_pedigree, selected_rows)        
        }
            
        nround <- nround + 1
        
        if(nround == maxn && prop_num > nrow(selected_df_by_pedigree)){
            new_row_remaining <- prop_num - nrow(selected_df_by_pedigree)
            take_all_pedigrees <- in_order[1:nrow(in_order), 1]
            sample_all_pedigree <- sample(take_all_pedigrees, new_row_remaining)
            new_val_data_remaining <- val_data_fr %>% anti_join(selected_df_by_pedigree, by = "ID")
            new_selected_df_remaining <- new_val_data_remaining %>% filter(Pedigree %in% sample_all_pedigree)
            selected_df_by_pedigree <- rbind(selected_df_by_pedigree, new_selected_df_remaining)
        }    
    }
    return(selected_df_by_pedigree)
 }   
