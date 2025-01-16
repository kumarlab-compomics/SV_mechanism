# TITLE OF PAPER

The objective is to describe structural (SVs) by features, including but not limited to local homology. 

This is done by:
- creating a workflow to describe local homology and
- training and applying an active unsupervised approach to describe underlying patterns of SVs.
- downstream analyses
  
The purpose of this is to better describe potential mechanisms of SV introduction, namely double strand break (DSB) repair. As such, SVs, their breakpoints and their flanking sequences are annotated with a number of features. These features include sequence, DNA accessibility, Replication Timing, repetitive elements, nonB DNA structures, and DNA shape. To robustly describe real SVs, we then generate 100 simulated SVs which match in the real SVs' length, type and chromosome. This therefore, allows for z-score standardization across all these features. 

## Required sources of data

### Structural Variants (SVs)

A vcf or csv containing resolved insertions and/or deletions of 50bp or greater. This file must contain a unique list of SVs with the following columns

-  CHROM; Chromosome
-  POS; Breakpoint position
-  REF and ALT; Reference and Alternative sequence
-  SVlen; SV length
-  SV_Type; SV type (accepted: INS, DEL, insertion, deletion)

### Human genome reference

Initial version was built using HG 38. 

### DNA accessibility

We annotated deleted regions, along with the flanks of insertions and deletions. Bigwigs which map to Chip-seq, DNase-seq and WGB-seq were gathered from ENCODE. For germline SVs, we used experiments conducted in H1-hESCs. The following experiments and accessibility markers were used

- ENCBS111ENC (DNase-seq, WGB-seq+, WGB-seq-); https://www.encodeproject.org/biosamples/ENCBS111ENC/ 
- ENCBS718AAA (CTCF, H2AFZ, H3K27ac, H3K27me3, H3K36me3, H3K4me1, H3K79me2, H3K9ac, H3K9me3, H4Kme1); https://www.encodeproject.org/biosamples/ENCBS718AAA/

### non-B DNA structures

To annotate non-B DNA structures in flanking regions, files from https://nonb-abcc.ncifcrf.gov/apps/ftp/browse were further processed by collapsing regions into bed files. 

### Repeat Masker motifs

In a similar fashion, the file https://www.repeatmasker.org/species/hg.html was processed into bed files to characterize known sites of LINE, Low complexity, LTRs, Satellites, simple repeats and SINEs. 

### Replication Timing 

Results from a 16-stage Repli-Seq experiment was downloaded from 4D nucleome _. This file was further processed using Repliseq to describe genomic locations by their S50. 

## SV annotation

Due to the size and scale of SVs in whole genome population-level studies, SVs from the vcf are split by chromosome, then further batched. Using HPC resources, each batch is then annotated automatically and in parallel using the execute_Master1.sh. From this master script, the following batch scripts, and python scripts, are triggered :

- execute_seqFeaturesSV.sh
  - adding_seqFeaturesSV_240118.py
- seqFeatures;
- flankSeq ;
- SVcoords ;
- RepliSeq ;
- searchRepeatMasker ;
- ENCODE ;

When running the master batch script, the following arguments must be provided : 

- Directory which holds the Project
- Subdirectly which holds the vcf
- Name of the vcf
- Directory to point to cell line used for DNA accessibility annotation

An example of the use of this command is: 
```
sbatch execute_Master1_240502.sh \
allHGSVC3_20240415_Freeze4_GRCh38 \
sv \
freeze4_GRCh38_sv_insdel_alt.vcf.SVlength.csv \
H1
```

## Simulation generation & Z-score calculation




