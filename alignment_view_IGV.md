# Starting

Indexing fist the reference. 

Please define which reference you will use:

a complete dataframe or a model from V9 region instead.

```bash
srun bowtie2-build ./reference/euk_V9.unique.fa ./refence/euk_V9unique
```

Then, lets align libs vs the indexed reference

```bash
for i in A2 A3 A4 A5 A10
do
srun bowtie2  -x ./reference/euk_V9.unique -q -1 04-X04-${i}-18S-AMB_*_L001_R1_001.fastq.gz -2 04-X04-${i}-18S-AMB_*_L001_R2_001.fastq.gz -S bwt_${i}.sam
done 

```

And compress into a binary format

```BASH
for i in A2 A3 A4 A5 A10
do
samtools view -bS bwt${i}.sam > bwt${i}.bam 
done 

```

And then sort its binary

```bash
for i in A2 A3 A4 A5 A10
do
samtools sort bwt${i}.bam -o bwt_${i}.bam.sorted
done
```

Let's index:

```bash
for i in A2 A3 A4 A5 A10
do
samtools index bwt_${i}.bam.sorted
done
```



Let's use `cufflinks` to make coordinates from the samples:

```bash
for i in A2 A3 A4 A5 A10
do
cufflinks --overlap-radius 1 \
             --library-type fr-firststrand \
             -o cufflinks.${i}.dir bwt_${i}.bam.sorted
done
```



and then merge each transcripts.gtf into matrix using cuffmerge:

> export PATH="$PATH:/home/rgomez/bin/cufflinks-2.2.1.Linux_x86_64/"

```bash
ls cufflinks.A*/transcripts.gtf > cufflinks.txt
cuffmerge -s reference/euk_V9.unique.fa cufflinks.txt

```



now copy to local machine and Use IGV:

```bash
scp rgomez@omica:/LUSTRE/bioinformatica_data/genomica_funcional/mothur/bowtie/*.sorted .
scp -r rgomez@omica:/LUSTRE/bioinformatica_data/genomica_funcional/mothur/bowtie/cufflinks.*.dir .
scp -r rgomez@omica:/LUSTRE/bioinformatica_data/genomica_funcional/mothur/bowtie/merged_asm .
scp rgomez@omica:/LUSTRE/bioinformatica_data/genomica_funcional/mothur/bowtie/reference/*.fa .

```





# visualization in IGV:
cufflinks.A10.dir
Open > File > Load FILE > Select within CUFFLINKS directory to load

Referece than use to align reads (genome/transcriptome assembly)
Merged *gtf file
Merge independent each *gtf file from condition (ex.cufflinks.ds.dir/Sp_ds.transcripts.gtf)
hits bam files from tophat (ex. Tophat.hs/hs_hits.bam)
Use your reference

/LUSTRE/bioinformatica_data/genomica_funcional/mothur/bowtie