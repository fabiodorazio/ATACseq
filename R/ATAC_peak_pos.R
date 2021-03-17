all.peaks.pgc <- ensembl.ids.atac(pgc.summit)
all.peaks.soma <- ensembl.ids.atac(soma.summit)

atac_list <- list(all.peaks.pgc, all.peaks.soma)

## get peaks coordinates for each promoter
.atac_positions <- function(y){
  max_duplicates <- ddply(y,.(geneId), nrow)
  max_duplicates <- max(max_duplicates$V1)
  
  # create columns to populate
  num_col_original <- ncol(y) + 1
  new_names <- paste0('Peak_n_',2:max_duplicates)
  y[,c(as.numeric(num_col_original):as.numeric(ncol(y)+max_duplicates-1))] <- 0
  colnames(y)[num_col_original:ncol(y)] <- new_names
  # ascending order
  y <- y[order(y$geneId),]
  # sort in ascending order within groups
  y <- y %>%
    group_by(geneId) %>%
    arrange(geneId, distanceFromProm) %>%
    data.frame()
  
  # fill the df with 6 empty rows to avoid data loss
  fake_df <- head(y)
  fake_df[,1] <- "f_ENSEMBL"
  fake_df[,2:ncol(fake_df)] <- 0
  y <- rbind(fake_df, y)
  
  # fill the data frame
  for(element in 1:nrow(y)){
    tryCatch({
      if(y$geneId[element] == y$geneId[element - 1] & y$geneId[element] != y$geneId[element - 2]){
        y$Peak_n_2[element-1] <- y$distanceFromProm[element]
      }
      else if(y$geneId[element] == y$geneId[element - 2] & y$geneId[element] != y$geneId[element - 3]){
        y$Peak_n_3[element-2] <- y$distanceFromProm[element]
      }
      else if(y$geneId[element] == y$geneId[element - 3] & y$geneId[element] != y$geneId[element - 4]){
        y$Peak_n_4[element-3] <- y$distanceFromProm[element]
      }
      else if(y$geneId[element] == y$geneId[element - 4] & y$geneId[element] != y$geneId[element - 5]){
        y$Peak_n_4[element-4] <- y$distanceFromProm[element]
      }
      else if(y$geneId[element] == y$geneId[element - 5] & y$geneId[element] != y$geneId[element - 6]){
        y$Peak_n_4[element-5] <- y$distanceFromProm[element]
      }
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
  }
  # remove top rows
  y <- y[-c(1:6),]
  # count repeated Peaks
  y3 <- y %>% 
    group_by(geneId) %>%
    dplyr::summarise(N = n()) %>%
    data.frame()
  
  # remove duplicated rows
  z <- y[!duplicated(y$geneId),]
  
  # add count to original table
  z$num_ATAC_peaks <- y3$N
  # reconsitute original table
  z$insideFeature <- NULL
  colnames(z)[2] <- 'Peak_n_1'
  
  # average distance
  z1 <- z
  z1[z1 == 0] <- NA
  z1$Avg_Peak_dist = apply(z1,1,function(x) mean(as.numeric(x[2:7])[!is.na(x[2:7])]))
  
  # add to original
  z$Avg_Peak_dist <- z1$Avg_Peak_dist
  
  return(z)
  
  #cage_obj_tata <- merge(cage_obj, z, by.x = "cluster", by.y = "sequence", all = TRUE)
}

atac_list_positions <- lapply(atac_list, .atac_positions)
