# TITLE OF PAPER

The objective is to describe structural (SVs) by features, including but not limited to local homology. 

This is done by:
- creating a workflow to describe local homology and
- training and applying an active unsupervised approach to describe underlying patterns of SVs.
- downstream analyses
  
The purpose of this is to better describe potential mechanisms of SV introduction, namely double strand break (DSB) repair. As such, SVs, their breakpoints and their flanking sequences are annotated with a number of features. These features include sequence, DNA accessibility, Replication Timing, repetitive elements, nonB DNA structures, and DNA shape. To robustly describe real SVs, we then generate 100 simulated SVs which match in the real SVs' length, type and chromosome. This therefore, allows for z-score standardization across all these features. 

## Step 0 : Real structural variant (SV) source

Start with a vcf with SVs. For example, you can download and place the integrated callset from HGSVC2 here:

```
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v2.0/integrated_callset/variants_freeze4_sv_insdel_alt.vcf.gz -P ./data/HGSVC2/
```

Next, this vcf should then be manipulated to contain: resolved insertions and/or deletions of 50bp or greater. This file must contain a unique list of SVs with the following columns:

-  CHROM; Chromosome
-  POS; Breakpoint position
-  REF and ALT; Reference and Alternative sequence
-  SVlen; SV length
-  SV_Type; SV type (accepted: INS, DEL, insertion, deletion)

As per our example, the file should be saved here:

```
./preprocess/HGSVC2/SVlen/sv/variants_freeze4_sv_insdel_alt.vcf.SVlength.csv
```

## Step 1 : Generating simulated SVs

Our goal is to use underlying repair features to better discriminate between DSB repair mechanisms. Since we are using SVs (50bp+), we noted that when calculating features such as mean GC content, the local or SV's GC content could depend on the SV length, chromosome of origin and the SV type. We therefore wanted to standardize numeric features. To do this rigorously, we generated 100 simulations for each real and unique SV. These simulations were generated using SURVIVOR. The simulations matched the real SVs in their: chromosome, length and type. 

The requirements for this command includes: 
- The "real" vcf/ csv which will be used to generate simulated SVs
- The SURVIVOR tool
- HG38 reference genome
  - This file should be stored as such : ```./data/Genome/hg38.fa ```
- Scratch space, if using an HPC.

To start the ``` execute_svSimsMaster.sh ``` script the following arguments must be provided:
1. The project name
2. sv (given we are generating simulated SVs held in the sv directory)
3. The name of the vcf-like file we created in Step 0. 

An example of the use of this command is from inside the ```./simulations``` directory should be : 

```
cd simulations

sbatch execute_svSimsMaster.sh \
HGSVC2 \
sv \
variants_freeze4_sv_insdel_alt.vcf.SVlength.csv
```

In the end you should have two large vcfs/csvs that represent the real and simulated SVs. The simulated csv should be 100x the size of the real csv. Therefore, as per the HGSVC file, you should have the simulated vcf-like file here:

```
./preprocess/HGSVC2/SVlen/svSIM/variants_freeze4_sv_insdel_alt.vcf.SVlength.Simulations.csv
```

## Step 2 : Annotating real and simulated SVs

Next, we want to annotate the real and simulated vcf-like files. [See more details about how and which features are annotated, along with data requirements below.]

Due to the scale and size of the whole-genome population-level studies, the vcf-like files are split by chromosome, then further batched into 10K SVs. For each chromosome, a set of bash scripts are deployed in parallel, which then starts bash and python scripts to annotate each batch. 

To start the ``` execute_Master1.sh ``` script the following arguments must be provided:
1. The project name
2. sv or svSIM, depending on the file we are annotating
3. The name of the file to annotate
4. The name of directory used to store epigenetic data 
   
As per our HGSVC2 example, the following should be deployed to annotate the real and simulated SVs. In this case, we will use H1-hESCs as the epigenetic proxy.

```
cd ../preprocess

sbatch execute_Master1.sh \
HGSVC2 \
sv \
variants_freeze4_sv_insdel_alt.vcf.SVlength.csv \
H1

sbatch execute_Master1.sh \
HGSVC2 \
svSIM \
variants_freeze4_sv_insdel_alt.vcf.SVlength.Simulations.csv \
H1
```

In ``` execute_Master1.sh ``` the following bash scripts are copied and customized to start the correct number of array jobs, as per the number of batches. Within these bash scripts, python scripts are then executed. 

- ``` execute_ENCODETemplate.sh``` : To calculate SV and flanking epigenetic profiles
  - ``` adding_epiFeaturesSV.py ``` 
  - ``` adding_epiFeaturesflanks.py ``` 
- ``` execute_RepliSeqTemplate.sh ``` : To annotate the replication timing of the SV's position
  - ``` adding_RepliSeq.py ```
- ``` execute_SVcoordsTemplate.sh ``` : This generates bed files to represent the up and downstream flanks to the SVs. We then annotate these beds based on their confidence to produce R-loops, Repeat Masker motifs and nonB DNA structures.
  - ``` adding_SVcoords.py ```
  - ``` adding_flankRLoop.py ```
  - ``` adding_flankRepeatMasker.py ```
  - ``` adding_flanknonBDNA.py ```
- ``` execute_chromoBandTemplate.sh ``` : To annotate the giemsa stain at an SV's position
  - ``` adding_chromoBand.py ```
- ``` execute_flankSeqTemplate.sh ``` : To extract reference genome sequences from the flanking coordinates. Next, we use these sequences to calculate sequence features, local homology (as per Blast alignments) and DNA shape
  - ``` adding_flankSeq.py ```
  - ``` adding_seqFeatures.py ```
  - ``` adding_BlastDNAShape.py ```
  - ``` adding_Blastmerges.py ```
- ``` execute_searchRepeatMaskerTemplate.sh ``` : To run Repeat Masker on the deleted and inserted sequences
  - ``` adding_searchRepeatMasker.py ```
  - ``` adding_searchRepeatMaskermerges.py ```
- ``` execute_seqFeaturesSVTemplate.sh ``` : To calculate sequence features of the SVs' sequences
  - ``` adding_seqFeaturesSV.py ```

Subject to successful annotation of each SV batch, we then must merge these annotations. We use the ``` execute_Master1Merge.sh ``` script, which requires the following arguments : 
1. The project name
2. sv or svSIM, depending on the file we are merging
3. The name of directory used to store epigenetic data 

As per our HGSVC2 example, here are the commands for this merging script for both the real and simulated SVs

```
sbatch execute_Master1Merge.sh \
HGSVC2 \
sv \
H1

sbatch execute_Master1Merge.sh \
HGSVC2 \
svSIM \
H1
```

### Data required for annotations

#### ``` execute_ENCODETemplate.sh ```

We annotated deleted regions, along with the flanks of insertions and deletions. We do so by using the ``` pyBigWig ``` module, along with a series of bigwigs provided by ENCODE. For our germline SVs, we use the H1-hESC cell line as a proxy for the epigenetic status of the genome. We begin by creating a directory here ``` ./data/ENCODE/H1 ```. We found in ENCODE, the Chip-Seq experiment and DNase-seq/WGB-seq experiment which were held H1 data was ENCBS718AAA and ENCBS111ENC, respectively. As such create the following directories

```
mkdir /data/ENCODE/H1/ENCBS111ENC 
mkdir /data/ENCODE/H1/ENCBS718AAA

cd ENCBS111ENC
mkdir DNase-seq WGB-Seq
cd ../ENCBS718AAA
mkdir CTCF H2AFZ H3K36me3 H3K9ac H3K9me3 H4K20me1 H3K4me1 H3K79me2 H3K27ac H3K27me3
```

For these experiments refer to the following locations to download the appropriate big wigs. Use wget to download the bigwigs in the appropriate directory

- For ENCBS111ENC : https://www.encodeproject.org/biosamples/ENCBS111ENC/
- For ENCBS718AAA : https://www.encodeproject.org/biosamples/ENCBS718AAA/

#### ``` execute_RepliSeqTemplate.sh ```

[ re-do description + instructions ] 

Results from a 16-stage Repli-Seq experiment was downloaded from 4D nucleome _. This file was further processed using Repliseq to describe genomic locations by their S50. 

We have provided the H1 processed Repli-Seq file here : ``` ./data/RepliSeq/H1/ ```


#### ``` execute_SVcoordsTemplate.sh ```

[ re-do description + instructions ] 

 - To annotate non-B DNA structures in flanking regions, files from https://nonb-abcc.ncifcrf.gov/apps/ftp/browse were further processed by collapsing regions into bed files.
 - In a similar fashion, the file https://www.repeatmasker.org/species/hg.html was processed into bed files to characterize known sites of LINE, Low complexity, LTRs, Satellites, simple repeats and SINEs. 

#### ``` execute_chromoBandTemplate.sh ```


#### ``` execute_flankSeqTemplate.sh ```

- basically all you need is the reference genome, which you already have


## Step 3 : Calculating z-scores across features

Now, that both the real and simulated SVs have been annotated, we will now match the real SVs to their matching simulations to calculate z-scores across the numeric features. We will use the ``` execute_Master1Zscore.sh ``` script to do this. The requirements are : 
1. The project name
2. sv
3. svSIM

Therefore, for the HGSVC2 example, we apply the following : 

```
sbatch execute_Master1Zscore.sh \
HGSVC2 \
sv \
svSim
```

## Step 4 : Applying the homology-based label to SVs

Upon generating the dataframes which hold features' z-scores, we can now begin labeling SVs. We start by applying strict thresholds based on the up-downstream flanking alignments. We consider both the length and 'quality' of the gap-less alignment we obtained from Blast. To do this, we use the ``` execute_mechIDsvSim.sh ``` script in the ``` ./analysis ``` directory. The requirements are as follows:

1. The project name
2. sv

For the HGSVC2 dataset, we would use the following : 

```
../analysis

sbatch execute_mechIDsvSim.sh \
HGSVC2 \
sv
```

We are left with the following classes of SVs: high local homology (HLH), intermediate local homology (ILH), no local homology (NLH) and 'Undefined'. 

Note: In this case, we are annotating the raw Blast coverage and identity scores, therefore, the bash script could be manipulated to use the pre-zscore generated files. 


## Step 5 : Identifying optimal clustering parameters, and training reusable models  




