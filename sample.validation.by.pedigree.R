sample.validation.by.pedigree <- function(val_pheno, val_pedigree_data, prop_num){
    
    val_data_fr <- val_pedigree_data %>% inner_join(val_pheno, by = "ID")
    in_order <- val_data_fr %>% count(Pedigree) %>% arrange(desc(n))
    more_than_one <- in_order %>% filter(n > 1)
    just_one <- in_order %>% filter(n == 1)

    if(prop_num <= nrow(more_than_one)){  # priority for pedigrees that have multiple lines
        selected_ped <- in_order[1:prop_num, 1]
        selected_df <- data.frame()
        for(i in selected_ped){    
            selected_rows <- sample_n(subset(val_data_fr, Pedigree == i), 1)    
            selected_df <- rbind(selected_df, selected_rows)
        }
        selected_df_out <- selected_df
    }else if(prop_num > nrow(more_than_one) && prop_num < nrow(in_order)){
        n_for_more_than_one <- nrow(more_than_one) # # number of peds with multiple lines
        n_for_just_one <- prop_num - n_for_more_than_one  # # number of peds with single lines
        just_one_ped_list <- just_one[1:nrow(just_one), 1]  # list of pedigrees with just one line
        just_one_ped_sample <- sample(just_one_ped_list, n_for_just_one)  # samples pedigree name to select

        selected_df_for_just_one <- val_data_fr %>% filter(Pedigree %in% just_one_ped_sample)

        selected_ped <- in_order[1:n_for_more_than_one, 1]
        selected_df_more_than_one <- data.frame()
        for(i in selected_ped){    
            selected_rows <- sample_n(subset(val_data_fr, Pedigree == i), 1)    
            selected_df_more_than_one <- rbind(selected_df_more_than_one, selected_rows)
        }

        selection_df_out <- rbind(selected_df_more_than_one, selected_df_for_just_one)

    }else if(prop_num > nrow(in_order)){  # select from peds with multiple lines multiple times

        n_for_more_than_one <- nrow(more_than_one) # number of peds with multiple lines
        selected_ped <- in_order[1:n_for_more_than_one, 1]
        selected_df_more_than_one <- data.frame()
        for(i in selected_ped){    
            selected_rows <- sample_n(subset(val_data_fr, Pedigree == i), 1)    
            selected_df_more_than_one <- rbind(selected_df_more_than_one, selected_rows)
        }

        n_for_just_one <- nrow(just_one) # number of peds with single lines
        just_one_ped_list <- just_one[1:n_for_just_one, 1]  # list of pedigrees with just one line
        just_one_df <- val_data_fr %>% filter(Pedigree %in% just_one_ped_list)

        first_dt <- rbind(selected_df_more_than_one, just_one_df)
        n_for_first_dt <- nrow(first_dt)
        n_for_second_round <- prop_num - n_for_first_dt

        second_round_df <- val_data_fr %>% anti_join(first_dt, by = "ID") # so no line is sampled twice

        selected_ped_second_round <- in_order[1:n_for_second_round, 1]
        selected_df_second_round <- data.frame()
        for(i in selected_ped_second_round){    
            selected_rows <- sample_n(subset(second_round_df, Pedigree == i), 1)    
            selected_df_second_round <- rbind(selected_df_second_round, selected_rows)
        }
        selected_df_out <- rbind(first_dt, selected_df_second_round)
    }
}   
