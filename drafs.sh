#Get beta diversity of plant survey

richness <- colSums(curated_table)
obs_beta <- nrow(curated_table)/mean(richness)

ggplot(method_statistics, aes(x=Curated,weights=Beta,fill=Curated)) +
  geom_bar(position="dodge") + 
  geom_hline(yintercept = obs_beta,linetype = 2) +
  facet_grid(Method ~Level) +
  xlab("") +
  ylab("Betadiversity (total richness / avg. plot richness)") +
  scale_fill_brewer(palette = "Set1") +
  theme_bw()+ theme(text = element_text(size=8)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

**Figure 6. Betadiversity of each method.**  
Betadiversity (calculated as total number of OTUs divided by the mean number of OTUs in the 130 sites) is plotted on the y-axis. Red bars represent the betadiversity of the taxonomic assignment from the un-curated methods and blue bars represent that of the curated methods. The dashed line indicates the betadiversity of the plant data (17.23) observed in the study for comparison, calculated in the same way. The “99/100%” clustering level denotes the pure DADA2 approach (100%) and SWARM with a d-value of 3 (99%). All methods initially have a higher betadiversity than can be expected from observational plant data – even CROP, which contain very few OTU. The betadiversity is reduced to realistic levels with curation for all methods.  

\pagebreak

 
```{r dpi=300,echo=FALSE}
tab_name <- file.path(main_path,"Table_OTU_match_rates.txt")
pident_frame <- read.table(tab_name, sep="\t", header=TRUE, as.is=TRUE)
pident_frame$level_pident <- 
 factor(pident_frame$level_pident,levels = c("99/100", "98.5", "98", "97", 
                                             "96", "95"))
pident_frame$method_pident <- 
 factor(pident_frame$method_pident,levels = c("VSEARCH", "SWARM", "DADA2(+VS)", 
                                              "CROP"))
#Violin plot of distribution of best matches, all methods
ggplot(pident_frame, aes(x=retained_or_discarded, 
                         y=pident,fill=retained_or_discarded)) +
  geom_violin() +
  facet_grid(method_pident~level_pident)+
  xlab("number of OTUs") +
  ylab("Best match on GenBank") +
  scale_fill_brewer(palette = "Set1") +
  theme_bw()+
  theme(text = element_text(size=8)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))


# MONDAY

/LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/data/run14_20190228_COI

sbatch dada2_coi.sh TESTNAME DATASET_PATH

# 
#!/bin/bash
#SBATCH --job-name=rdp
#SBATCH -N 2
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24

echo "Fecha inicio: `date`"
echo "Ejecutandose con $SLURM_JOB_CPUS_PER_NODE"
echo "Numero de nodos: $SLURM_NNODES y CPU por nodo: $SLURM_CPUS_ON_NODE"
echo "CPUs totales= $SLURM_NPROCS"
echo "Los nodos utilizados son: $SLURM_NODELIST"

# to run:
# sbatch rdp_assign.sh run014_20190304_COI_t1/run014_t1_ASVs.fasta 99 outdir

# Variables

fasta=$1
boots=$2
outdir=$3

DB=/LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/dbs/
DB_REF="midori_unique_DB-0.fasta"
DB_TAX="midori_unique_DB-0.tax"

mothur="/LUSTRE/bioinformatica_data/genomica_funcional/bin/mothur/mothur"

echo "Asignacion taxonomia de $fasta con $DB_REF and $DB_TAX"

cd $SLURM_SUBMIT_DIR

echo "Starting ... !"


$mothur "#system(mkdir -p $outdir);set.dir(output=$outdir, tempdefault=$DB);summary.seqs(fasta=$fasta, processors=$SLURM_NPROCS);classify.seqs(fasta=current, reference=$DB_REF, taxonomy=$DB_TAX, iters=1000, cutoff=$boots);get.current();quit()"

exit

#

/LUSTRE/bioinformatica_data/genomica_funcional/MG_COI/dada2_asv/run14_20190228_COI/dada2_asv


# 1. Sample_ID,
# 2. Sample_Name,
# 3. Sample_Plate,
# 4. Sample_Well,
# 5. I7_Index_ID,
# 6. index,
# 7. I5_Index_ID,
# 8. index2,
# 9. Sample_Project,
# 10. Description

# Sample_ID code: 
# <###>_<CRUICE>_<STATION>_<MARKER>_<GROUP>
# Ej. 001_X4_A1_18S_AMB

file=run006_SampleSheet.csv 
paste -d, <(awk 'NR>21' $file| cut -d',' -f1) \
          <(awk 'NR>21' $file| cut -d',' -f1 | sed 's/_/,/g') \
          <(awk 'NR>21' $file | cut -d',' -f6,8) > metadata.csv

#
for file in *csv;
do
paste -d, <(awk 'NR>21' $file| cut -d',' -f1) \
          <(awk 'NR>21' $file| cut -d',' -f1 | sed 's/_/,/g') \
          <(awk 'NR>21' $file | cut -d',' -f6,8)
done > metadata.csv


Ophiura cf. ljungmani


