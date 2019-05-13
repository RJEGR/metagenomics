# cigom_coi_midori2worms.R
# Febrero 2019
# Miguel Martinez | marmigues@gmail.com
#
# R functions to convert midori taxonomy to WoRMS taxonomy framework
# names2worms() - Translates a names vector to WoRMS taxonomy

## Names conversion ----
# Modified from "cigom_coi_fish_convert.R > names2worms"

"names2worms" <- function(mynames) {
    # 'mynames' should be a vector with unique names (no NA's)
    # output is a DF with the columns in 'fields_interest'

    # Require
    require(stringr)
    require(worrms)
    require(taxize)
    
    # Oficial taxonomy ranks
    TL <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    
    # Taxonomy levels as DF
    TL_df <- data.frame(rank=factor(TL, levels =TL, ordered = T) )
    
    # Fields in WoRMS-Record that are of interest 
    fields_interest <- data.frame(fields=c("AphiaID", "scientificname", "status", "rank", 
                                           "valid_AphiaID", "valid_name", "valid_authority", 
                                           "isMarine", "isBrackish", "isFreshwater", "isTerrestrial", "modified"), stringsAsFactors = F) 
    
    # Columns for converted-DF
    info_names <- c("ori_name", "AphiaID", "rank", "status", "modified",
                    paste0(TL, "_wm"), paste0(TL, "_wid"), "valid_AphiaID", "valid_name", "valid_authority",
                    "isMarine", "isFreshwater", "isTerrestrial")
                    #"isBrackish", queda fuera

    # Object for converted names
    n_info <- data.frame(matrix(ncol=25))
    # Loop IDs to retrieve info
    for(n in 1:length(mynames)){    
        ## Worms
        # Get ID (taxize)
        out_id <- get_wormsid(mynames[n], accepted = F, ask=T)
        # If only Species names appear (or lower levels), can use one, script will correct the name
        # Evaluate ID existence
        if(is.na(out_id) | is.null(out_id)){
            # AphiaID doesnt exist
            n_info[n,] <- c(mynames[n], rep(NA, 24))
        } else {
            # AphiaID exists, now check Record
            # Get record using ID
            out_rec <- wm_record(as.integer(out_id))
            # Transform to DF
            out_rec_df <- data.frame(do.call(rbind, out_rec), stringsAsFactors = F)
            # Get info of interest
            fields_df <- merge(x=fields_interest, y=out_rec_df, by.x="fields", by.y=0, all.x=T)
            names(fields_df)[2] <- "to_validate"
            # Evaluate Record
            # If there's no valid-AphiaID or the Status is 'quarantined' or 'deleted'
            if(is.na(fields_df[which(fields_df$fields %in% "valid_AphiaID"), "to_validate"]) | fields_df[which(fields_df$fields %in% "status"), "to_validate" ] %in% c("quarantined", "deleted")){
                # AphiaID cant be helped
                # Generate a line filled with NAs
                n_info[n, ] <- c(mynames[n], rep(NA, 24))
            } else {
                # AphiaID is ok, a trusted record can be obtained using the valid-AphiaID
                # retrieve good-record from good-ID
                good_aphia <- as.integer(fields_df[which(fields_df$fields %in% "valid_AphiaID"), "to_validate"])
                sec_rec <- wm_record(good_aphia)
                # Transform to DF
                sec_rec_df <- data.frame(do.call(rbind, sec_rec), stringsAsFactors = F)
                # Get info of interest
                good_sec_rec <- merge(x=fields_df, y=sec_rec_df, by.x="fields", by.y=0, all.x=T)
                names(good_sec_rec)[3] <- "updated"
                
                # Retrieve classification from good-ID
                sec_class <- data.frame(wm_classification(good_aphia), stringsAsFactors = F)
                # Evaluate duplicated ranks
                if(nrow(sec_class) != length(unique(sec_class$rank))){
                    # Keep only first rank ocurrence 
                    sec_class[duplicated(sec_class$rank), ] <- NA
                }#/ if dup-ranks
                # Combine with ranks of interest
                sec_class <- merge(x=TL_df, y=sec_class, by="rank", all.x = T)
                sec_class <- sec_class[order(sec_class$rank, decreasing = F), ]
        
                # Fill DF with recovered info of interest
                n_info[n, ] <- c(
                    mynames[n], 
                    # info from Record (4 fields)
                    good_sec_rec[which(good_sec_rec$fields %in% "AphiaID"), "updated"],
                    good_sec_rec[which(good_sec_rec$fields %in% "rank"), "updated"],
                    good_sec_rec[which(good_sec_rec$fields %in% "status"), "updated"],
                    # Keep only date of modification (remove time)
                    sub("T.*", "", good_sec_rec[which(good_sec_rec$fields %in% "modified"), "updated"]),                    
                    # info from Taxonomy (7 names and 7 AphiaIDs)
                    sec_class[which(sec_class$rank %in% c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")), "scientificname"],
                    sec_class[which(sec_class$rank %in% c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")), "AphiaID"],
                    # more info from Record (6 fields)
                    good_sec_rec[which(good_sec_rec$fields %in% "valid_AphiaID"), "updated"],
                    good_sec_rec[which(good_sec_rec$fields %in% "valid_name"), "updated"],
                    good_sec_rec[which(good_sec_rec$fields %in% "valid_authority"), "updated"],
                    good_sec_rec[which(good_sec_rec$fields %in% "isMarine"), "updated"],
                    good_sec_rec[which(good_sec_rec$fields %in% "isFreshwater"), "updated"],
                    good_sec_rec[which(good_sec_rec$fields %in% "isTerrestrial"), "updated"])
            }#/ ifelse  eval-record  
        }#/ ifelse eval-ID
    }#/ for-Names

    # Add names
    names(n_info) <- info_names
    
    # Return WORMS information about mynames
    assign("n_info", n_info, .GlobalEnv)
    
} #fin-names2worms()

"names2worms_" <- function(mynames, accepted = TRUE, ask = FALSE) {
    # 'mynames' should be a vector with unique names (no NA's)
    # output is a DF with the columns in 'fields_interest'

    # Require
    require(stringr)
    require(worrms)
    require(taxize)
    
    # Oficial taxonomy ranks
    TL <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    
    # Taxonomy levels as DF
    TL_df <- data.frame(rank=factor(TL, levels =TL, ordered = T) )
    
    # Fields in WoRMS-Record that are of interest 
    fields_interest <- data.frame(fields=c("AphiaID", "scientificname", "status", "rank", 
                                           "valid_AphiaID", "valid_name", "valid_authority", 
                                           "isMarine", "isBrackish", "isFreshwater", "isTerrestrial", "modified"), stringsAsFactors = F) 
    
    # Columns for converted-DF
    info_names <- c("ori_name", "AphiaID", "rank", "status", "modified",
                    paste0(TL, "_wm"), paste0(TL, "_wid"), "valid_AphiaID", "valid_name", "valid_authority",
                    "isMarine", "isFreshwater", "isTerrestrial")
                    #"isBrackish", queda fuera

    # Object for converted names
    n_info <- data.frame(matrix(ncol=25))
    # Loop IDs to retrieve info
    for(n in 1:length(mynames)){    
        ## Worms
        # Get ID (taxize)
        # out_id <- get_wormsid(mynames[n], accepted = F, ask=T)
        out_id <- get_wormsid(mynames[n], accepted = accepted, ask = ask)
        # If only Species names appear (or lower levels), can use one, script will correct the name
        # Evaluate ID existence
        if(is.na(out_id) | is.null(out_id)){
            # AphiaID doesnt exist
            n_info[n,] <- c(mynames[n], rep(NA, 24))
        } else {
            # AphiaID exists, now check Record
            # Get record using ID
            out_rec <- wm_record(as.integer(out_id))
            # Transform to DF
            out_rec_df <- data.frame(do.call(rbind, out_rec), stringsAsFactors = F)
            # Get info of interest
            fields_df <- merge(x=fields_interest, y=out_rec_df, by.x="fields", by.y=0, all.x=T)
            names(fields_df)[2] <- "to_validate"
            # Evaluate Record
            # If there's no valid-AphiaID or the Status is 'quarantined' or 'deleted'
            if(is.na(fields_df[which(fields_df$fields %in% "valid_AphiaID"), "to_validate"]) | fields_df[which(fields_df$fields %in% "status"), "to_validate" ] %in% c("quarantined", "deleted")){
                # AphiaID cant be helped
                # Generate a line filled with NAs
                n_info[n, ] <- c(mynames[n], rep(NA, 24))
            } else {
                # AphiaID is ok, a trusted record can be obtained using the valid-AphiaID
                # retrieve good-record from good-ID
                good_aphia <- as.integer(fields_df[which(fields_df$fields %in% "valid_AphiaID"), "to_validate"])
                sec_rec <- wm_record(good_aphia)
                # Transform to DF
                sec_rec_df <- data.frame(do.call(rbind, sec_rec), stringsAsFactors = F)
                # Get info of interest
                good_sec_rec <- merge(x=fields_df, y=sec_rec_df, by.x="fields", by.y=0, all.x=T)
                names(good_sec_rec)[3] <- "updated"
                
                # Retrieve classification from good-ID
                sec_class <- data.frame(wm_classification(good_aphia), stringsAsFactors = F)
                # Evaluate duplicated ranks
                if(nrow(sec_class) != length(unique(sec_class$rank))){
                    # Keep only first rank ocurrence 
                    sec_class[duplicated(sec_class$rank), ] <- NA
                }#/ if dup-ranks
                # Combine with ranks of interest
                sec_class <- merge(x=TL_df, y=sec_class, by="rank", all.x = T)
                sec_class <- sec_class[order(sec_class$rank, decreasing = F), ]
        
                # Fill DF with recovered info of interest
                n_info[n, ] <- c(
                    mynames[n], 
                    # info from Record (4 fields)
                    good_sec_rec[which(good_sec_rec$fields %in% "AphiaID"), "updated"],
                    good_sec_rec[which(good_sec_rec$fields %in% "rank"), "updated"],
                    good_sec_rec[which(good_sec_rec$fields %in% "status"), "updated"],
                    # Keep only date of modification (remove time)
                    sub("T.*", "", good_sec_rec[which(good_sec_rec$fields %in% "modified"), "updated"]),                    
                    # info from Taxonomy (7 names and 7 AphiaIDs)
                    sec_class[which(sec_class$rank %in% c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")), "scientificname"],
                    sec_class[which(sec_class$rank %in% c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")), "AphiaID"],
                    # more info from Record (6 fields)
                    good_sec_rec[which(good_sec_rec$fields %in% "valid_AphiaID"), "updated"],
                    good_sec_rec[which(good_sec_rec$fields %in% "valid_name"), "updated"],
                    good_sec_rec[which(good_sec_rec$fields %in% "valid_authority"), "updated"],
                    good_sec_rec[which(good_sec_rec$fields %in% "isMarine"), "updated"],
                    good_sec_rec[which(good_sec_rec$fields %in% "isFreshwater"), "updated"],
                    good_sec_rec[which(good_sec_rec$fields %in% "isTerrestrial"), "updated"])
            }#/ ifelse  eval-record  
        }#/ ifelse eval-ID
    }#/ for-Names

    # Add names
    names(n_info) <- info_names
    
    # Return WORMS information about mynames
    assign("n_info", n_info, .GlobalEnv)
    
} #fin-names2worms_()

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
      mc_df <- names2worms_(mc_names)
      
      # Successful online conversion
      mc_df_ok <- mc_df[!is.na(mc_df$AphiaID), ]
      
      # outList <- list( "mc_df"  = mc_df)
      # return(outList)
      
      # }}} # aqui se puede acabar, pero continuamos ...
      # Still-non-converted:
      
      snc_names <- mc_df[is.na(mc_df$AphiaID), "ori_name"]
      
      outList <- list( "mc_df_ok" = mc_df_ok, "mc_df"  = mc_df, "snc_names" = snc_names )
      
      return(outList)
      
    } } } #fin-ranks2worms()

