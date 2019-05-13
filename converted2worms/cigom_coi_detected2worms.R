# cigom_coi_detected2worms.R
# Febrero 2019
# Miguel Martinez | marmigues@gmail.com
#

# how to: Rscript --vanilla cigom_coi_detected2worms.R wang.taxonomy midori2worms_conversion.tab

# Script that converts taxonomy from mothur-rdp-midori assigned ASVs to WoRMS nomenclature.
# It keeps record of detected/converted names in a `conversion.tab` file.

# wang.taxonomy file is DB_tax specific
# midori_unique_DB_0 taxonomy file has 2 columns (tab-sep, no header):
# ASV - char - ID
# taxonomy - char - Taxonomy string (semicolon-sep) with 8 levels: root, superkingdom, K, P, C, O, F, G, S

# NOTAS:
# 1. What to do with ?
    # S Homo sapiens 
    # S Canis lupus
    # S Camponotus floridanus
    # S Liposcelis entomophila
    # F Dermestidae
    # 
    # Create a black list to remove them? prepare fake profiles?
    # ver nc_lineages

# Load libraries ----
.cran_packages <- c("stringr", "tidyverse")
.ropensci_packages <- c("worrms", "taxize")

# 1. cran_package
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst], dep=TRUE, repos='http://cran.us.r-project.org')
}

# 2. ropensci_package
.inst <- .ropensci_packages %in% installed.packages()

if (any(!.inst)) {   
  if (!require('devtools', lib.loc = lib.loc)) {
    install.packages("devtools", dep=TRUE, repos='http://cran.us.r-project.org')
  } else 
    sapply(c(.ropensci_packages), install_github, username = 'ropensci')
}

# 3. jamba
if (!require("jamba")) {   
  if (!require('devtools', lib.loc = lib.loc)) {
    install.packages("devtools", dep=TRUE, repos='http://cran.us.r-project.org')
  } else 
    devtools::install_github("jmw86069/jamba")
}

# Load packages into session, and print package version
sapply(c(.ropensci_packages, .cran_packages, 'ggplot2', 'jamba'), require, character.only = TRUE)
# sapply(c(.ropensci_packages, .cran_packages), packageVersion)


## Prepare space/variables ----

## Clean workspace
rm(list=ls()); 

## commit 2: args insertion<< ----
# # # 
# Args
# # #

args = commandArgs(trailingOnly=TRUE)


if (length(args)<2) {
  stop("!!!\n. 
       PLEASE, INPUT BOTH FILES IN THE SYNTAXIS AS FOLLOW EXAMPLE:\n
       Rscript --vanilla cigom_coi_detected2worms.R wang.taxonomy midori2worms_conversion.tab .\n", call.=FALSE) 
  } else { 
    wangtax = args[1]
    conv_file = args[2] }


path <- '/Users/cigom/metagenomics/COI/run014'
setwd(path)
## Input: wang.taxonomy
wangtax = 'run014_t2_ASVs.midori_unique_DB_0.wang.taxonomy' # 'run012_relax_ASVs.midori_unique_DB_0.wang.taxonomy_99'
# Input: Conversion file (Names already converted)
conv_file <- "~/Documents/GitHub/metagenomics/converted2worms/midori2worms_conversion.tab"


## Variables
# Oficial taxonomy ranks
TL <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# colnames in conversion file (conv_file)
conv_names <- c("ori_name", "AphiaID", "rank", "status", "modified", 
                paste0(TL, "_wm"), paste0(TL, "_wid"), "valid_AphiaID", "valid_name", "valid_authority", 
                "isMarine", "isFreshwater", "isTerrestrial")
                #"isBrackish", queda fuera
#
# OUTPUT 
#
## commit 3: outcomes defined << ----

path <- getwd()

out_dir <- paste0(path, '/asv_analisis')

tag <- strsplit(wangtax, "[.]")[[1]]

# make out_dir
system(command = paste0("mkdir -p ", out_dir), intern = F)

#Out_name
out_tax <- paste0(out_dir, '/', tag[1], '_',tag[2], '.',tag[3],'.wm.taxonomy')

## Non-converted lineages
nc_tax <- paste0(out_dir,'/', tag[1],'_nonconverted_lineages.txt')

# 
# Functions 
## Load home-made functions
#
# To load file:  wang.taxonomy
source(file = "~/Documents/GitHub/metagenomics/converted2worms/cigom_load_dada2_tables.R")
# To convert names
source(file = "~/Documents/GitHub/metagenomics/converted2worms/cigom_coi_midori2worms.R")


## Make reference backup copy ----

system(command = paste0("cp ", conv_file, " ", conv_file, ".", format(Sys.time(), "%Y%b%d")))    
#^-  modificar Sys.time por un comando que obtenga la fecha de ultima modificacion o de creacion

## Load already converted ----

# Load
conv_df <- read.table(conv_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)
dim(conv_df)
# 16292     25  # v1 (20 Febrero 2019)

## Load assignation result ----

# wang.taxonomy, threshold = 99
rdp_df <- load.mothurdp_boots(wangtax, 80)

# 13,872    8 #run012
# 21,278   8 # run014

# Remove unclassified and replace underscores by spaces
rdp_df <- data.frame(rdp_df) %>% 
            mutate_all(funs(str_replace_all(., c(".+_unclassified"=NA, "_"=" "))))

rdp_df$ASV <- str_replace(rdp_df$ASV, " ", "_")

# Tag starting level (SL) (aka calculate resolution)
rdp_df$SL <- 7
for(t in 7:2){
    rr <- t-1  #rank real 
    rc <- t+1  #column-index (rank) para corregir
    if(t == 7){
        rdp_df$SL[is.na(rdp_df[,rc])] <- rr
    } else {
        rdp_df$SL[is.na(rdp_df[,rc]) & rdp_df$SL <= t ] <- rr
    }
}
# Add column for converted level (CL)

table(rdp_df$SL)
#    2    3    4    5    6    7 (run12)
# 8666   96  214  308  396 4192 

# 2     3     4     5     6     7  (run14)
# 12432    98   521   211   467  7549 

## Loop levels for conversion ----
## commit 4: determine the undetermined << ---- 
##  DF to trim converted rows

n_undetermined <-  nrow(rdp_df[rdp_df$Phylum == 'Undetermined',])

cat('\n Removing', n_undetermined,  'Undetermined Phylum from the dataset', '\n')

if (n_undetermined > 0) { rdp_new <- rdp_df[-which(rdp_df$Phylum %in% "Undetermined"), ] 
} else {rdp_new = rdp_df} 

dim(rdp_new)

# DF to add conversions (all Undetermined)
rdp_ready <- rdp_df[which(rdp_df$Phylum %in% "Undetermined"), ] #Modify this DF by adding converted rows
#rdp_ready$CL <- 0 #no aun
dim(rdp_ready)
# 8000 9

# non-converted names
nc_lineages <- character()

# rdp_new <- rdp_new[sample(nrow(rdp_new), 50), ] # to test

x <- apply(rdp_new[TL[-1]], 2, ranks2worms, conv_df = conv_df)

mc_df <- do.call(rbind, lapply(x, `[[`, 'mc_df'))
mc_df_ok <- do.call(rbind, lapply(x, `[[`, 'mc_df_ok'))
ac_df <- do.call(rbind, lapply(x, `[[`, 'ac_df'))
snc_names <- sapply(x, `[[`, 'snc_names')

#.mc_df_ = do.call(rbind, sapply(lapply(x, `[[`, 'mc_df'), nrow))
.mc_df_ = data.frame(sapply(lapply(x, `[[`, 'mc_df'), nrow))
.mc_df_ok_ = data.frame(sapply(lapply(x, `[[`, 'mc_df_ok'), nrow))
.ac_df_ = data.frame(sapply(lapply(x, `[[`, 'ac_df'), nrow))
.scn_names_ = data.frame(sapply(sapply(x, `[[`, 'snc_names'), length))

y <- data.frame(matrix(0, nrow = length(TL[-1]), ncol = 4 ))
rownames(y) <- TL[-1]

y[which(rownames(.mc_df_) %in% rownames(y)),1] <- .mc_df_
y[which(rownames(.mc_df_ok_) %in% rownames(y)),2] <- .mc_df_ok_
y[which(rownames(.ac_df_) %in% rownames(y)),3] <- .ac_df_
y[which(rownames(.scn_names_) %in% rownames(y)),4] <- .scn_names_

names(y) <- c('mc_df', 'mc_df_ok', 'ac_df', 'snc_names')


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


#L <- 2 # Level testing
# Loop TL levels
for(L in 7:2){
    cat(paste0("\nAnalyzing level: ", TL[L], "\n"))
    
    # Identify unique names
    u_names <- sort(na.omit(unique(rdp_new[,TL[L]])))
    length(u_names)
    # 494   L7
    # 112   L6
    # 32    L5
    # 16    L4
    # 9     L3
    # 4     L2
    
    # Verificacion que asegure que al menos hay un nombre unico para convertir
    if(length(u_names) > 0){
        
        # Obtain already-converted names (ac)
        ac_df <- conv_df[which(conv_df$ori_name %in% u_names), ]
        dim(ac_df)
        # 487 25    L7
        # 0   25    L6
        # 0   25    L5
        # 0   25    L4
        # 0   25    L3
        # 0   25    L2
    
        # Missing conversion (mc)
        mc_names <- sort(setdiff(u_names, ac_df$ori_name))
        length(mc_names)
        # 7     L7
        # 112   L6
        # 32    L5
        # 16    L4
        # 9     L3
        # 4     L2
    
        # Verificacion que asegure que al menos hay un nombre para convertir
        if(length(mc_names) > 0 ){
          
## Commit 5: input logical vars in names2worms to non-interactive name selection << ----               
          
            # Online conversion (using homemade function)
            mc_df <- names2worms_(mc_names) # To non-interactive name use names2worms_ function
            # Successful online conversion
            mc_df_ok <- mc_df[!is.na(mc_df$AphiaID), ]      
            dim(mc_df_ok)
            # 0  25     L7
            # 108 25    L6
            # 27    25  L5
            # 14    25  L4
            # 9     25  L3
            # 4     25  L2
                # Add this result to already converted (conv_df)
                conv_df <- rbind(conv_df, mc_df_ok)                                    # Aqui se incrementa la DB de nombres ya convertidos, imprimir a archivo al terminar el loop
        
                dim(conv_df)
                # 16,293  25    L7
                # 16,401  25    L6
                # 16,428  25    L5
                # 16,442  25    L4
                # 16,451  25    L3
                # 16,455  25    L2

            # Still-non-converted
            snc_names <- mc_df[is.na(mc_df$AphiaID), "ori_name"]
            length(snc_names)
            # 7     L7
            # 4     L6
            # 5     L5
            # 2     L4
            # 0     L3
            # 0     L2
        
                # Keep track of non converted lineages (esta linea esta mal debido a que el objecto nc_df solo recupera la ultima lista del ciclo)
                nc_df <- rdp_new[which(rdp_new[,TL[L]] %in% snc_names), TL] %>% unite("lineage", TL, sep=";", remove = TRUE)
                nc_lineages <- c(nc_lineages, unique(nc_df$lineage))

        }#/if hay nombres para names2worms
    
        # Merge successful conversions 'ac' + 'mc'
        conv_ok <- rbind(ac_df, mc_df_ok)
        dim(conv_ok)
        # 487 25    L7
        # 108 25    L6
        # 27  25    L5
        # 14  25    L4
        # 9   25    L3
        # 4   25    L2

        # Integrate successful conversions (SC) to starting DF
        SC_df <- merge(x=rdp_new, y=conv_ok[, c("ori_name", paste0(TL, "_wm"))], by.x=L+1, by.y="ori_name", all.x = TRUE)
        # clean DF
            #SC_df <- SC_df[!is.na(SC_df[,paste0(TL[L-2], "_wm")]), ]                   # No estoy seguro de que columna verificar para eliminar lo no encontrado
            SC_df <- SC_df[!is.na(SC_df$Kingdom_wm), ]
        # format
        SC_df <- SC_df[, c("ASV", paste0(TL, "_wm"), "SL")]
        names(SC_df) <- c("ASV", TL, "SL")
        dim(SC_df)
        # 4163  9   L7
        # 419   9   L6
        # 307   9   L5
        # 219   9   L4
        # 100   9   L3
        # 666   9   L2
        
        # Add to 'ready' df
        rdp_ready <- rbind(rdp_ready, SC_df)
        dim(rdp_ready)
        # 12,163 9   L7
        # 12,582 9   L6
        # 12,889 9   L5
        # 13,106 9   L4
        # 13,206 9   L3
        # 13,872 9   L2
    
        table(rdp_ready$SL)
        #    2    3    4    5    6    7 
        # 8000    0    0  305  396 4188   #L5
        # 8000    0  214  308  396 4188   #L4
        # 8000   96  214  308  396 4192   #L3
        # 8666   96  214  308  396 4192   #L2
        
        # Identify level-above in still-non-converted ori_names
        rdp_new <- rdp_new[-which(rdp_new$ASV %in% SC_df$ASV), ]
        dim(rdp_new)    
        # 1,709   9  L7
        # 1,290   9  L6
        # 983     9  L5
        # 766     9  L4
        # 666     9  L3
        # 0       9  L2
    
        table(rdp_new$SL)
        #   2   3   4   5   6   7 
        # 666  96 214   3   0   4   #L5
        # 666  96   0   0   0   4   #L4
        # 666   0   0   0   0   0   #L3
        # --                        #L2
    
    }#/ if hay nombres unicos
    
} #/ for, looping TL




# TODO FUNCIONA HASTA AQUI :::: COMO HACER EL MERGE RANK BY RANK ???


# Integrate successful conversions (SC) to starting DF
# SC_df <- merge(x=rdp_new, y=conv_ok[, c("ori_name", paste0(TL, "_wm"))], by.x=L+1, by.y="ori_name", all.x = TRUE)
# SC_df <- merge(x=rdp_new, y=conv_ok[, c("ori_name", paste0(TL, "_wm"))], by.x=TL[-1], by.y="ori_name", all.x = TRUE)

# clean DF
#SC_df <- SC_df[!is.na(SC_df$Kingdom_wm), ]

# format
#SC_df <- SC_df[, c("ASV", paste0(TL, "_wm"), "SL")]
#names(SC_df) <- c("ASV", TL, "SL")
#dim(SC_df)

# Add to 'ready' df
#rdp_ready <- rbind(rdp_ready, SC_df)
# continuar con el flujo del ciclo !!!!
dim(rdp_ready)

# Two important DFs are:
dim(rdp_ready)  # Converted taxfile-DF
# 13,872    9
dim(conv_df)    # Increased/Updated reference DF
# 16,455    9

## Get something to plot ---

cp_ready <- rdp_ready[-which(rdp_ready$Phylum %in% "Undetermined"),]
num_undet <- nrow(rdp_ready) - nrow(cp_ready)

# Tag converted level (CL) (aka calculate resolution)
cp_ready$CL <- 7

for(t in 7:2){
    rr <- t-1  #rank real 
    rc <- t+1  #column-index (rank) para corregir
    if(t == 7){
        cp_ready$CL[is.na(cp_ready[,rc])] <- rr
    } else {
        cp_ready$CL[is.na(cp_ready[,rc]) & cp_ready$CL <= t ] <- rr
    }
}

table(cp_ready$CL)
#    2    3    4    5    6    7 
#  666  101  209  261  290 4345 

## Plots ----

## Cual es la distribucion de la resolucion taxomica?
L_tables <- data.frame(Rank=TL[2:7], Level=2:7, SL=table(cp_ready$SL), CL=table(cp_ready$CL))
# Wide to long
L_tables_g <- gather(L_tables[,c(1,2,4,6)], key="Stage", value = "ASVs", -c("Rank", "Level"))
L_tables_g$Stage <- factor(L_tables_g$Stage, levels=c("SL.Freq", "CL.Freq"), labels=c("Initial", "Converted"), ordered = T)
L_tables_g$Rank <- factor(L_tables_g$Rank, levels = TL[2:7], ordered = T)

# Lineplot. Distribucion de las asignaciones Inicial y final

ggplot(L_tables_g, aes(x=Rank, y=ASVs, group=Stage, fill=Stage, color=Stage)) + 
    #geom_col(position=position_dodge()) +
    geom_point(size=2, alpha=0.6) + geom_line(size=1, alpha=0.6) +
    scale_color_brewer(palette = "Set1") +
    labs(x="", title="Midori to WoRMS taxonomy conversion", 
               subtitle = paste0("Taxonomy resolution distribution.\nTotal ASVs = ", nrow(cp_ready), " ('Undetermined' removed = ", num_undet, ")"),
               caption = paste0("Formating ", tag[1], " file"))

## Que nivel taxonomia tienen los ASVs al terminar mejor o peor?
cp_ready$L_diff <- cp_ready$CL - cp_ready$SL
# L_diff = 0 Se convirtio el mismo nivel
# L_diff > 0 El nivel convertido es menor al inicial (ganancias)
# L_diff < 0 El nivel convertido es mayor al inicial (perdidas)
L_diff_g <- data.frame(table(cp_ready$L_diff))
# -2   -1    0    1    2    4 
# 40  160 5641   25    2    4
# Barplot
ggplot(L_diff_g, aes(x=Var1, y=Freq)) + 
    geom_col(width = 0.4, fill="thistle4")+ #"indianred") + 
    geom_text(aes(label=Freq), size=3, vjust=-0.2) +
    labs(x="Levels changed", y="ASVs (n)", 
         title="Midori to WoRMS taxonomy conversion", 
         subtitle=paste0("Won (+) and lost (-) levels\nTotal ASVs = ", nrow(cp_ready)),
         caption = paste0("Formating ", tag[1], " file")) +
    theme()

## Clean and print converted result ----    

# Keep good columns
wm_ready <- rdp_ready[,c("ASV", TL)]
# Paste taxonomy string
wm_ready <- wm_ready %>% unite("taxonomy", TL, sep=";", remove = T)
# Este paso se puede evitar si los ASV tuvieran la misma cantidad de digitos
wm_ready <- mixedSortDF(wm_ready, byCols = 1, decreasing = FALSE) #{jamba}
# Print out
write.table(wm_ready, file = out_tax, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
## Save Increased/Updated reference DF ----

write.table(conv_df, conv_file, sep = "\t", row.names = FALSE, col.names = TRUE)

## Print conflicting lineages ----
## Commit 6: comment if conflicts in nc_linege << ----

if(length(nc_lineages)>0){
    write.table(nc_lineages, file=nc_tax, sep="\t", append = TRUE, row.names = FALSE, col.names = FALSE)
} else {
  cat('There were not lineage converted!')
  }

## Workspace Backup ----

# Pick where left
#rm(list=ls()); load(file="/home/fatallis/cigom/scripts/mx_pipeline/asv_analisis/cigom_coi_detected2worms.RData")

## Commit 7: Save RData in the asv_analisis folder  << ----

# Guardar workspace
save.image(file = paste0(out_dir,'/', tag[1],'_detected2worms.RData'))
#save.image(file= "/home/fatallis/cigom/scripts/mx_pipeline/asv_analisis/cigom_coi_detected2worms.RData")

cat('\nConversion done!\n')

quit(save = 'no')
