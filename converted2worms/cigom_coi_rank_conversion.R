
# apply(rdp_newTL[, TL], 2, ranks2worms)
# length_minus <- function (x) {length(x) - 5000}
# Ej. apply(rdp_new[, TL], 2, length_minus)
# ranks2worms(rdp_new$Phylum, conv_df)
# x <- apply(rdp_new[TL[-1]], 2, ranks2worms, conv_df = conv_df, id.names = rdp_new$ASV)
# rownames(rdp_new) <- rdp_new$ASV
# test_rdp_new <- head(rdp_new, n = 100)
# # x <- apply(test_rdp_new[TL[-1]], 2, ranks2worms, conv_df = conv_df)

"ranks2worms" <- function(x, conv_df) {
# x = Undetermined Phylum removed from the dataset
# conv_df = conversion file

  #rdp_new_ <- data.frame(rank = rdp_new$Phylum, ASV = id.names, stringsAsFactors=FALSE)
  rdp_new <- x
  # Identify unique names
  u_names <- sort(na.omit(unique(rdp_new)))
  
  # Verificacion que asegure que al menos hay un nombre unico para convertir
  if(length(u_names) > 0){   
      
      # Obtain already-converted names (ac)
      ac_df <- conv_df[which(conv_df$ori_name %in% u_names), ]
      
      # Missing conversion (mc)
      mc_names <- sort(setdiff(u_names, ac_df$ori_name))
      
      # Verificacion que asegure que al menos hay un nombre para convertir
      if(length(mc_names) > 0 ){

      # Online conversion (using homemade function)
      mc_df <- names2worms(mc_names)
     
       # Successful online conversion
      mc_df_ok <- mc_df[!is.na(mc_df$AphiaID), ]

     # outList <- list( "mc_df"  = mc_df)
     # return(outList)
      
     # }}} # aqui se puede acabar, pero continuamos ...
     # Still-non-converted:
      
      snc_names <- mc_df[is.na(mc_df$AphiaID), "ori_name"]
      
      outList <- list( "mc_df_ok" = mc_df_ok, "mc_df"  = mc_df, "snc_names" = snc_names, "ac_df" = ac_df )
      # es necesario incluir el archivo
      return(outList)

      } } }

# rdp_new <- rdp_new[sample(nrow(rdp_new), 100), ] # to test

x <- apply(rdp_new[TL[-1]], 2, ranks2worms, conv_df = conv_df)

# mc_df_ok <- lapply(x, `[[`, 'mc_df_ok')
mc_df <- do.call(rbind, lapply(x, `[[`, 'mc_df'))
mc_df_ok <- do.call(rbind, lapply(x, `[[`, 'mc_df_ok'))
ac_df <- do.call(rbind, lapply(x, `[[`, 'ac_df'))
snc_names <- sapply(x, `[[`, 'snc_names')

# Proportion of convertion:

# null_zero <- function(x) {x[which(sapply(x, nrow) ==  'NULL')] <- 0  }

y <- data.frame(matrix(0, nrow = length(TL[-1]), ncol = 4 ))
rownames(y) <- TL[-1]

.mc_df_ = do.call(rbind, sapply(lapply(x, `[[`, 'mc_df'), nrow))
.mc_df_ok_ = do.call(rbind, sapply(lapply(x, `[[`, 'mc_df_ok'), nrow))
.ac_df_ = do.call(rbind, sapply(lapply(x, `[[`, 'ac_df'), nrow))
.scn_names_ = data.frame(sapply(sapply(x, `[[`, 'snc_names'), length))

y[which(rownames(.mc_df_) %in% rownames(y)),1] <- .mc_df_
y[which(rownames(.mc_df_ok_) %in% rownames(y)),2] <- .mc_df_ok_
y[which(rownames(.ac_df_) %in% rownames(y)),3] <- .ac_df_
y[which(rownames(.scn_names_) %in% rownames(y)),4] <- .scn_names_

names(y) <- c('mc_df', 'mc_df_ok', 'ac_df', 'snc_names')
y

for(L in 7:2){ nc_df <- list(rdp_new[which(rdp_new[, TL[L]] %in% snc_names[[TL[L]]]), TL]) }

# nc_df <- do.call(rbind, lapply(nc_df, data.frame,  stringsAsFactors=FALSE)) %>% unite("lineage", TL, sep=";", remove = TRUE)


nc_df <- rbind(
  rdp_new[which(rdp_new[, TL[7]] %in% snc_names[[TL[7]]]), TL],
  rdp_new[which(rdp_new[, TL[6]] %in% snc_names[[TL[6]]]), TL],
  rdp_new[which(rdp_new[, TL[5]] %in% snc_names[[TL[5]]]), TL],
  rdp_new[which(rdp_new[, TL[4]] %in% snc_names[[TL[4]]]), TL],
  rdp_new[which(rdp_new[, TL[3]] %in% snc_names[[TL[3]]]), TL],
  rdp_new[which(rdp_new[, TL[2]] %in% snc_names[[TL[2]]]), TL]) 



#  keep only unique/distinct lineage.

nc_df %>% 
  unite("lineage", TL, sep=";", remove = TRUE) %>% 
  distinct() ->  nc_lineages 


# Se incrementa la DB de nombres ya convertidos:
conv_df <- rbind(conv_df, mc_df_ok)

# Merge successful conversions 'ac' + 'mc'
conv_ok <- rbind(ac_df, mc_df_ok)

dim(conv_ok)

# TODO FUNCIONA HASTA AQUI :::: COMO HACER EL MERGE RANK BY RANK ???

# Integrate successful conversions (SC) to starting DF
# SC_df <- merge(x=rdp_new, y=conv_ok[, c("ori_name", paste0(TL, "_wm"))], by.x=L+1, by.y="ori_name", all.x = TRUE)

SC_df <- merge(x=rdp_new, y=conv_ok[, c("ori_name", paste0(TL, "_wm"))], by.x=TL[-1], by.y="ori_name", all.x = TRUE)

# clean DF
SC_df <- SC_df[!is.na(SC_df$Kingdom_wm), ]

# format
SC_df <- SC_df[, c("ASV", paste0(TL, "_wm"), "SL")]
names(SC_df) <- c("ASV", TL, "SL")
dim(SC_df)

# Add to 'ready' df
rdp_ready <- rbind(rdp_ready, SC_df)
dim(rdp_ready)


table(rdp_ready$SL)

# Identify level-above in still-non-converted ori_names
rdp_new <- rdp_new[-which(rdp_new$ASV %in% SC_df$ASV), ]
dim(rdp_new)    


table(rdp_new$SL)


}#/ if hay nombres unicos

  