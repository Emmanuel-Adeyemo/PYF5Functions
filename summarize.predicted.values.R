summarize.predicted.values <- function(out_list){
    
    twt_values <- out_list[[2]]
    twt_values <- twt_values %>% arrange(ID)
    rownames(twt_values) <- twt_values[,1]
    twt_values <- twt_values[,-1]
    twt_av_values <- rowMeans(twt_values, na.rm = T)
    twt_av_values <- as.data.frame(twt_av_values)
    colnames(twt_av_values) <- "pred_twt"
    
    vsk_values <- out_list[[3]]
    vsk_values <- vsk_values %>% arrange(ID)
    rownames(vsk_values) <- vsk_values[,1]
    vsk_values <- vsk_values[,-1]
    vsk_av_values <- rowMeans(vsk_values, na.rm = T)
    vsk_av_values <- as.data.frame(vsk_av_values)
    colnames(vsk_av_values) <- "pred_vsk"
    
    dis_values <- out_list[[4]]
    dis_values <- dis_values %>% arrange(ID)
    rownames(dis_values) <- dis_values[,1]
    dis_values <- dis_values[,-1]
    dis_av_values <- rowMeans(dis_values, na.rm = T)
    dis_av_values <- as.data.frame(dis_av_values)
    colnames(dis_av_values) <- "pred_dis"
    
    out <- cbind(twt_av_values, vsk_av_values, dis_av_values)
    
    return(out)    
}
