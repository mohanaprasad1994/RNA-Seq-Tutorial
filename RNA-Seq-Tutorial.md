# RNA-Seq Pipeline

The overall pipieline is the following:
![Pipeline](pipeline.png|alt=octocat)

We start with the RNA sequence read files in .bam format format. We use HTSeq to count gene expressions for each sample. We then compile the results for all samples. If your samples come from multiple populations, you might want to analyse each of the population separately for EQTL. We use DESeq2 to variance stabilize each of the populations. We then normalize the variance stabilized counts to zero mean and unit variance. We then run Surrogate Variable Analysis to find hiiden covariates that need to be corrected for. Finally we run EQTL for each population with their corresponding genotypes, normalized-variance stabilized gene counts, SVs and other known covaiates to be corrected for.

Each of the steps are explained in detail below.

# HTSeq

HTSeq takes the raw read files (.bam files) and the .gtf file and gives the expression of each gene in a sample.

Check out http://htseq.readthedocs.io/en/master/count.html and set the -s, -m, -i,-t and -r fields correctly. 

#### Scripts for running HTSeq

##### To run HTSeq for a single sample: 
run_htseq.sh: script to run htseq-count with the required parameters
```sh
$ ./run_htseq.sh <.bam filepath> <.gtf filepath> <output filepath>
```
>run_htseq.sh  
```sh
#!/bin/bash

module load python/2.7
module load HTSeq/0.9.1

htseq-count -f bam  -r pos -t gene -i gene_id -m union -s reverse \
    $1  \
    $2  > \
    $3
```
##### To run htseq on entire dataset:
Set the Global variables in run_htseq.py appropriately and run on scg4:
```sh
$ module load python/2.7
$ python run_htseq.py
```
It takes about 3-4 hrs per sample
>run_htseq.py
```sh
import os
import glob

# path to the .gtf or .gff file 
GTF_FILE_PATH = "/srv/gsfs0/projects/montgomery/mohanas/africa/rna/ref/gencode.v27.annotation.chr-renamed.gtf"

# path to folder containing all reads. Assumption: the .bam files are located like "<ID>/Aligned.sortedByCoord.out.bam" relative to this folder. example ID = ERR188033
READS_DIR = "/srv/gsfs0/projects/montgomery/mdegorte/data/Africa/RNA/globus/SJ.nU2.nM0.maxO0.nSamples2.all.tab/"

# Path to the directory to store the results. Make sure it exists.
RESULTS_DIR = "/srv/gsfs0/projects/montgomery/mohanas/africa/rna/htseq-output/stranded/"

# list all the starting letters of the folders that need to be processed. 
START_IDs = ["HG", "GM", "NA"]

# path to script that htseq-count
SCRIPT_PATH = "/srv/gsfs0/projects/montgomery/mohanas/africa/rna/script/run-htseq/run_htseq.sh"

ct = 0
for folder in glob.glob(READS_DIR + "*"):
    Id = folder.strip().split('/')[-1]
    flag = False
    for name in START_IDs:
        if name in Id:
            flag = True
    if flag:
        print Id
        ct += 1
        read_file = READS_DIR + Id + "/Aligned.sortedByCoord.out.bam"
        out_file = RESULTS_DIR + Id + ".ReadsPerGene.out.tab"
        cmd = "qsub -cwd -M mohanas@stanford.edu -m ae -e EOs -o EOs -l h_vmem=8G " + SCRIPT_PATH + " " + read_file + " " + GTF_FILE_PATH + " " + out_file
        print cmd 
        os.system(cmd)
        #if ct >= 3:
        #    break

print ct
```

#### Merge all the HTSeq output into one matrix
To compile all htseq output files into one file:

Set the global variables appropriately and run:
```sh
$ module load python/2.7
$ python compile_gene_reads.py
```

>compile_gene_reads.py
```py
import os
import pandas as pd
import glob

GENE_READS_DIR = "/srv/gsfs0/projects/montgomery/mohanas/africa/rna/htseq-output/stranded/"
OUTPUT_FILE_PATH = "/srv/gsfs0/projects/montgomery/mohanas/africa/rna/GeneReadsMatrixStranded.csv"

empty_file_list = []
matrix = {}
ct = 0
for file_name in glob.glob(GENE_READS_DIR + "*"):
    print file_name
    Id = file_name.strip().split('.')[0].split('/')[-1]
    matrix[Id] = {}
    fil = open(file_name, 'r')
    for line in fil:
        key = line.strip().split()[0]
        val = line.strip().split()[1]
        matrix[Id][key] = val
    print Id, len(matrix[Id].keys())
    if len(matrix[Id].keys()) == 0:
        empty_file_list.append(Id)
    ct += 1
    #if ct == 3:
    #    break 

df = pd.DataFrame(matrix)
df.to_csv(OUTPUT_FILE_PATH)
print df
print "List of Ids with empty input files: ", empty_file_list
```


### Note:
- Rename samples here if need be 
- Make sure the chromosomes are names similarly in the read files and the .gtf files ( 1,2,3,... or Chr1, Chr2, Chr3, ….)
  To remove “chr” from .gtf file and rename "chr1" to "1":
    ```sh
    sed 's/^chr//' A.gtf > B.gtf
    ```


# Variance Stabilization with DESeq2
Remove genes with no expression (zeros in all individuals) before running variance stabilization.
To run for each population individually: Just do the same steps on the subset of the count
_data and the meta-data files

```R

## To install DESeq
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")

## Load csv
count_data <- read.csv("~/gene_expression/data/GeneReadsMatrixStrandedRenamed.csv", header=TRUE, row.names = 1 )
count_data = as.matrix(count_data)

## Remove unwanted rows
count_data <- count_data[ ! rownames(count_data) %in% c("__alignment_not_unique", "__ambiguous", "__no_feature", "__not_aligned", "__too_low_aQual"),]

## Remove genes that are not expressed in any of the samples
count_data = count_data[ rowSums(count_data)!=0, ] 

head(count_data)



## coldata <- read.csv("~/gene_expression/data/AGR_600RNA_sample_metadata_CP_300317.csv", header=TRUE, row.names=1)

## Meta data file with one entry per sample. row.names is the column number for index 
coldata <- read.csv("~/gene_expression/data/AGR_600RNA_sample_metadata_CP_300317_599Samples.csv", header=TRUE, row.names=1)

### For each population individually

# YRI_coldata = coldata[coldata$Population == 'YRI',]
# YRI_count_data = count_data[,rownames(YRI_coldata)]
# now perform the rest of the analysis on these variables.

## Shuffle so that the ordering of samples in coldata and count_data are the same
all(rownames(coldata) %in% colnames(count_data)) #True

all(rownames(coldata) == colnames(count_data)) # FALSE

count_data <- count_data[, rownames(coldata)]
all(rownames(coldata) == colnames(count_data)) # TRUE

## Perform Variance Stabilization
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = coldata, design = ~  1) 
dds=estimateSizeFactors(dds)
countData_lsc <- counts(dds, normalized=T)
normData=assay(varianceStabilizingTransformation(dds))

## Normalize the Variance Stabilized data
res = t(scale(t(normData)))

write.csv(res, file="~/gene_expression/data/deseq2-scaled-genecounts.csv")

```

# SVA

We can use SVA to correct for known effects as well
You need to give normalized variance stabilization output for SVA
If you want to do SVA over the entire dataset, then you can keep the variable of interest as "Population" (see ```mod``` in code). We can also use SVA to correct for known effects (```mod0```). We generally use "be" method to find the right number of SVs to be corrected for. You can also mention arbitary number of SVs to be corrected for.

If you are doing SVA for each population separately: 
1. Use the "two-step" method instead of  the default "irw"
2. This method is a little different and doesn't need ```mod0```. ```mod``` is just ```~1```
3. This method also determines the number of SVs by itself. It cannot give arbitary number of SVs




```R
#refer: http://bioconductor.org/packages/release/bioc/vignettes/sva/inst/doc/sva.pdf
# refer: http://jtleek.com/genstats/inst/doc/02_13_batch-effects.html
## To install sva
source("https://bioconductor.org/biocLite.R")
biocLite("sva")

library(sva)

### For population specific
# YRI_mod = model.matrix(~1, data=YRI_coldata)
# YRI_svobj = sva(YRI_res,mod=YRI_mod,method=c("two-step"))
# YRI_SVs <- cbind(c(row.names(YRI_coldata)), YRI_svobj$sv)
# write.csv(YRI_SVs, file = "~/gene_expression/data/YRI_surrogates_7.csv", row.names=FALSE)



mod = model.matrix(~as.factor(Population), data=coldata)
mod0 = model.matrix(~1,data=coldata)

#number of factors to be estimated
n.sv = num.sv(res,mod,method="be")

#Apply sva to find surrogate variables
svobj = sva(res,mod,mod0,n.sv=n.sv)

# Plot the probabilities
hist(svobj$pprob.gam, 
     main="Histogram for posterior probability that each gene is associated with one or more latent
variables", 
     xlab="Probability of association with surrogate variables", 
     ylab="Number of genes",
     border="blue", 
     col="green")

hist(svobj$pprob.b, 
     main="Histogram for posterior probability that each gene is associated with population", 
     xlab="Probability of association with population", 
     ylab="Number of genes",
     border="blue", 
     col="green")
#Adjust for surrogate variables

#without correcting for surrogate variables
pValues = f.pvalue(res,mod,mod0)
qValues = p.adjust(pValues,method="BH")

#Correction for surrogate variables
modSv = cbind(mod,svobj$sv)
mod0Sv = cbind(mod0,svobj$sv)
pValuesSv = f.pvalue(res,modSv,mod0Sv)
qValuesSv = p.adjust(pValuesSv,method="BH")

#print # of genes 
sum(qValuesSv < 0.05)
sum(qValues < 0.05)

# write surrogates to csv
SVs <- cbind(c(row.names(coldata)), svobj$sv)
write.csv(SVs, file = "~/gene_expression/data/surrogates.csv", row.names=FALSE, column.names=FALSE)

# to correct for SVs (if needed to be done)
RES=residuals(lm(res ~ SVs))
```

# EQTL Analysis:
EQTL is generally done for each population separately. Run variance stabilization, SVA for each of the population separately and then do EQTL on each of them.

**Input**: Normalized DESeq2 output ( normalized variance stabilized gene counts matrix),  SVs from SVA to be corrected for,  genotype information (vcf files)

**Output**: 

We use FastQTL for our EQTL analysis. Download the latest version of FastQTL from http://fastqtl.sourceforge.net/ 
It requires three files - the phenotypes file (.bed.gz + index files), the genotypes file (.vcf.gz + index files) and the covariates file (.txt).

We will now see how to generate each of them in detail.


#### Generating Phenotype file

**Input**: Normalized output of DESeq2 variance stabilized gene counts matrix, GTF file 
**Output**: Tab delimited .bed file with the following format
```
#Chr    start   end	ID	Sample_ID1  Sample_ID2  Sample_ID3  Sample_ID4
1	685395	685396	ENSG456 -1.13	1.18	-0.03	0.11
1	700304	700305	ENSG789 -1.18	1.32	-0.36	1.26
```

- Merge GTF file and your gene counts.
- Sort the file based on chromosomes and start index to get phenotypes.bed
    ```py
    GTF_FILE_PATH = "gencode.v27.annotation.chr-renamed.gtf"
    SCALED_DESEQ2_OUT_FILE_PATH = "deseq2-scaled-genecounts.csv" # Make sure the chromosomes are named similarly in the phenotypes and the vcf file
    OUTPUT_BED_FILE_PATH = "phenotypes.bed"
    
    # Read and strore the required parts of .gtf file
    gtf_dict = {}
    gtf_file = open(GTF_FILE_PATH, 'r')
    
    for i in range(5):
        gtf_file.readline()
        
    ct = 0
    chr_ids = []
    for line in gtf_file:
        l = line.strip().split()
        if l[2] == "gene": ### filter only genes
            gene_id = l[9][1:-2]
            chr_id = l[0]
            start = l[3]
            end = l[4]
            #print gene_id, chr_id, start, end
            gtf_dict[gene_id] = [chr_id, start, end, gene_id]
            chr_ids.append(chr_id)
            ct += 1
    
    print ct
    print set(chr_ids) # print all chromosome ids - just for checking
    
    # Load and display the deseq2 output
    scaled_deseq2_out_df = pd.read_csv(SCALED_DESEQ2_OUT_FILE_PATH, index_col=0)
    display(scaled_deseq2_out_df)
    
    # merge the .gtf file parts and the deseq2 output
    bed_file = open(OUTPUT_BED_FILE_PATH, 'w')
    header = "#" + "\t".join(["Chr", "start", "end", "ID"] + scaled_deseq2_out_df.columns.values.tolist())
    print header
    bed_file.write(header.strip() + "\n")
    
    ct = 0
    lines = []
    for row in scaled_deseq2_out_df.iterrows():
        ID, data = row
        #if gtf_dict[ID][0] not in ['X', 'Y', 'M']:
    
        lines.append(gtf_dict[ID] + data.tolist())
        ct += 1
    print ct
    lines = sorted(lines, key=lambda x: (str(x[0]),int(x[1]))) ### needs to be sorted for indexing to work
    
    for i,line in enumerate(lines):
        to_write = "\t".join(str(e) for e in line)
        if i == len(lines) - 1:
            bed_file.write(to_write.strip())
        else:
            bed_file.write(to_write.strip() + "\n")
    bed_file.close()
    
    ```
- Index the .bed file : ```bgzip phenotypes.bed && tabix -p bed phenotypes.bed.gz```


#### Generating Genotypes file:

- The ID column should have unique Ids. If not, do the following
-- To unzip .gz to a new location: ``` gunzip < file.vcf.gz > /new_location/file.vcf ```
-- To rename ID column and make it unique: ```bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' myfilename.vcf > new.vcf```
-- Compress and re-index: ``` bgzip new.vcf && tabix -p vcf new.vcf.gz```
-- Scripts that do the above: Run add_unique_ids.py after setting the global variables accordingly
    add_unique_ids.sh
    ```sh
    #!/bin/bash

    module load python/2.7
    module load bcftools/1.3.1
    
    gunzip < $1 > $2
    
    bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' $2 > $3
    
    rm -f $2
    
    /srv/gsfs0/projects/montgomery/mohanas/tools/tabix-0.2.6/bgzip $3 && /srv/gsfs0/projects/montgomery/mohanas/tools/tabix-0.2.6/tabix -p vcf $4

    rm -f $3
    ```
    add_unique_ids.py
    ```sh
    import os
    import glob
    import pandas
    
    READS_DIR = "/srv/gsfs0/projects/montgomery/mdegorte/data/Africa/RNA/globus/AGRVCFs/"
    SCRIPT_PATH = "/srv/gsfs0/projects/montgomery/mohanas/africa/rna/script/get_genotypes/add_unique_ids.sh"
    # Path to the directory to store the results. Make sure it exists.
    RESULTS_Dir = "/srv/gsfs0/projects/montgomery/mohanas/africa/genotypes/"
    ct = 0
    for file_path in glob.glob(READS_DIR + "*"):
        l = file_path.strip().split('/')[-1].split('.')
        if "no_mask" in l and "with_related" in l and l[-1] == 'gz' and l[-2] == 'vcf':
    #        print file_path
            ct += 1
            print str(ct) + "  " + file_path
            file_name = file_path.strip().split('/')[-1][:-3]
            param_1 = file_path
            param_2 = RESULTS_Dir + file_name
            param_3 = RESULTS_Dir + file_name[:-4] + ".id_filled.vcf"
            param_4 = param_3 + ".gz"
            cmd = "qsub -cwd -M mohanas@stanford.edu -m ae -e EOs -o EOs -l h_vmem=150G " + SCRIPT_PATH + " " + param_1 + " " + param_2 + " " + param_3 + " " + param_4
            print cmd
            os.system(cmd)
    ```
    It takes about 2-4 hrs per chromosome file
    
- At times, we might have to liftover genotypes files from different reference genome versions. Here we show an example of lifting over from hg19 to hg38. We will be using CrossMap (http://crossmap.sourceforge.net/ ) for this purpose. It can liftover many file formats including .vcf, .gtf, .bed, etc.
-- For lifting over .vcf files, you need input_chain_file in .chain or .chain.gz format and "ref_genome_file" (genome sequence file of 'target assembly') in FASTA foramt.
-- download the appropriate .chain and .fa files from http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/ and http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/ . You can find them from http://hgdownload.cse.ucsc.edu/downloads.html#liftover 
-- Command to run: ``` CrossMap.py vcf hg19ToHg38.over.chain.gz ../../genotypes/original_vcf/22.minAC1.no_mask.with_related.vcf.gz hg38.fa ../../genotypes/22.minAC1.no_mask.with_related.hg38.vcf > non_mkk.chr22.liftover.log 2>&1 ```
-- You might have to sort the vcf file before you run sequence number renaming, zipping and indexing. ``` cat X_PAR1.minAC1.no_mask.with_related.id_filled.hg38.vcf | vcf-sort > X_PAR1.minAC1.no_mask.with_related.id_filled.hg38.sorted.vcf ``` 


- or use picard (https://broadinstitute.github.io/picard/command-line-overview.html#LiftoverVcf ). You might want to create a sequence dictionary of the fasta file before calling liftover (https://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary )
-- ``` java -jar ../../../tools/picard.jar CreateSequenceDictionary R=hg38.fa O=hg38.fa.dict ```
-- ``` java -jar /srv/scratch/mohanas/tools/picard.jar LiftoverVcf I=/srv/scratch/mohanas/africa/genotypes/mkk_vcf/chr21.vcf.gz O=chr21.hg38.vcf CHAIN=/srv/scratch/mohanas/africa/rna/ref/hg19ToHg38.over.chain.gz REJECT=chr21.hg38.rejected_variants.vcf R=hg38.fa ```

- Merge genotypes files from multiple sources:

#### Generating Covariates file:

We would be using the output from SVA
Use the following code to get the covariates.txt file in the proper format
```py
SVs_FILE_PATH = "surrogates_100.csv"
NUM_SVs = 100

sv_df = pd.read_csv(SVs_FILE_PATH, index_col=0)
sv_df.columns = ["SV" + str(i+1) for i in range(NUM_SVs)]
sv_df.index.names = ['id']
display(sv_df.T)
sv_df.T.to_csv("covariates_"+str(NUM_SVs)+".txt", sep='\t', index_label="id")
```
Now zip the file: ``` bgzip covariates_100.txt ```

Note: You might also want to add known covariates like RIN, Gender, Pool Number for correction.

### To Run
```
/srv/persistent/bliu2/tools/fastqtl/bin/fastQTL --vcf /srv/scratch/mohanas/africa/genotypes/21.minAC1.no_mask.with_related.id_filled.vcf.gz --bed /srv/scratch/mohanas/africa/rna/phenotypes/YRI_phenotypes.bed.gz --out EQTL-Results/21.minAC1.no_mask.with_related.id_filled.YRI_phenotypes.YRI_covariates_7 --cov /srv/scratch/mohanas/africa/rna/covariates/YRI_covariates_7.txt.gz --threshold 0.05 --region 21 --exclude-samples /srv/scratch/mohanas/africa/rna/exclude_samples/YRI.exc
```

#### Note:
- VCF file - make sure the sample names in the vcf file match the sample names in phenotypes file.
- Make sure the chromosomes are named similarly in the phenotypes and the vcf file. (Either both use 1,2,3,.. Or Chr1, Chr2, ….)
- Make sure all the Genotype and Phenotype files have the corresponding index files
    To create index files:
    ```sh
    bgzip phenotypes.bed && tabix -p bed phenotypes.bed.gz
    bgzip genotype.vcf && tabix -p vcf genotype.vcf.gz
    ```
- FastQTL expects uniques values in each row for the ID column of genotypes file. (Sometimes it is all “.” in vcf files. Change it to say “Snap_<chr #>_<line #>”). You might have to unzip, change this column, zip and index all the vcf files again
- Sometimes, you might want to shrink the .vcf file by filtering only the samples needed in your analysis. Use the following sample command where ```sample_names.txt``` contains list of sample IDs, one per line. You have to zip and index after this step.
    ```
    vcftools --gzvcf /srv/gsfs0/projects/montgomery/mohanas/africa/genotypes/21.minAC1.no_mask.with_related.id_filled.vcf.gz --out /srv/gsfs0/projects/montgomery/mohanas/africa/genotypes/21.minAC1.no_mask.with_related.id_filled.599_samples --keep /srv/gsfs0/projects/montgomery/mohanas/africa/rna/sample_names.txt --recode
    ```
    


# Differential Gene Expression with DESeq2
As a separate analyse we do deseq2 on all pairs of population
1) You should filter genes to ignore psuedo-genes and keep only the required ones. (you can get gene types from the .gtf file)
2) Decide the threshold at which you are going to call a gene 'expressed', e.g. a gene is expressed if it has an average expression >5 counts, and zero counts in no more than 20% of your individuals (to avoid tails). This is an ad-hoc step and it depends on your sample size. 
3) Create new count matrix, only with expressed genes.

Sample code for Differential Expression Analysis:
```R
count_data <- read.csv("~/gene_expression/data/GeneReadsMatrixStrandedRenamed.csv", header=TRUE, row.names = 1 )
count_data = as.matrix(count_data)
count_data <- count_data[ ! rownames(count_data) %in% c("__alignment_not_unique", "__ambiguous", "__no_feature", "__not_aligned", "__too_low_aQual"),]
count_data = count_data[ rowSums(count_data)!=0, ] 
head(count_data)

coldata <- read.csv("~/gene_expression/data/AGR_600RNA_sample_metadata_CP_300317_599Samples.csv", header=TRUE, row.names=2)


all(rownames(coldata) %in% colnames(count_data)) #True
all(rownames(coldata) == colnames(count_data)) # FALSE
count_data <- count_data[, rownames(coldata)]
all(rownames(coldata) == colnames(count_data)) # TRUE

library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = count_data, colData = coldata, design = ~Gender + RIN + PoolNumber + Population) # by default the last variable in design is the condition
dds$Population <- relevel(dds$Population, ref = "MKK") # reference level is MKK
dds <- DESeq(dds)

resultsNames(dds)

res <- results(dds, contrast=c("Population","MSL","MKK"), alpha = 0.05, pAdjustMethod = "BH")
resLFC <- lfcShrink(dds, coef="Population_MSL_vs_LWK")
plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))

# find # genes with adjusted p-value less than a value
summary(res)
sum(res05$padj < 0.05, na.rm=TRUE)


# refer  https://shiring.github.io/rna-seq/deseq2/teaching/2016/09/29/DESeq2-course

summary(res)
mcols(res)$description

# order results table by the smallest adjusted p value:
res <- res[order(res$padj),]

library("dplyr")
library(ggplot2)
library(ggrepel)

results = as.data.frame(mutate(as.data.frame(res), sig=ifelse(res$padj<0.05, "FDR<0.05", "Not Sig")), row.names=rownames(res))
head(results)
DEgenes_DESeq <- results[which(abs(results$log2FoldChange) > log2(1.5) & results$padj < 0.05),]
p = ggplot2::ggplot(results, ggplot2::aes(log2FoldChange, -log10(pvalue))) +
  ggplot2::geom_point(ggplot2::aes(col = sig)) +
  ggplot2::scale_color_manual(values = c("red", "black")) +
  ggplot2::ggtitle("Volcano Plot of DESeq2 analysis")
p + ggrepel::geom_text_repel(data=results[1:10, ], ggplot2::aes(label=rownames(results[1:10, ])))

# plot counts
par(mfrow=c(1,3))
for (i in 1:3){
  gene <- rownames(res)[i]
  main = gene
  DESeq2::plotCounts(dds, gene=gene, intgroup="Population", main = main)
}

```

# Extras

Commands to Run:

```sh 
../../tools/FastQTL-2.165.linux/bin/fastQTL.1.165.linux --vcf /srv/gsfs0/projects/montgomery/mohanas/africa/genotypes/21.minAC1.no_mask.with_related.id_filled.vcf.gz --bed phenotypes_mini.bed.gz --region 21 --out EQTL-Results/21.minAC1.no_mask.with_related.id_filled.phenotypes_mini.covariates_34.txt.gz  --cov covariates_34.txt.gz --threshold 0.05
```

```sh 
qsub -cwd -M mohanas@stanford.edu -m ae -e EOs -o EOs -l h_vmem=200G   run_fastQTL.sh /srv/gsfs0/projects/montgomery/mohanas/africa/genotypes/21.minAC1.no_mask.with_related.id_filled.vcf.gz phenotypes_mini.bed.gz 21 EQTL-Results/21.minAC1.no_mask.with_related.id_filled.phenotypes_mini.covariates_34.txt.gz  /srv/gsfs0/projects/montgomery/mohanas/africa/rna/covariates_34.txt.gz 
```

To get subset of columns in vcf
need to zip and index after this
```
vcftools --gzvcf /srv/gsfs0/projects/montgomery/mohanas/africa/genotypes/21.minAC1.no_mask.with_related.id_filled.vcf.gz --out /srv/gsfs0/projects/montgomery/mohanas/africa/genotypes/21.minAC1.no_mask.with_related.id_filled.599_samples --keep /srv/gsfs0/projects/montgomery/mohanas/africa/rna/sample_names.txt --recode
```

```
qsub -cwd -M mohanas@stanford.edu -m ae -e EOs -o EOs -l h_vmem=200G /srv/gsfs0/projects/montgomery/mohanas/africa/rna/script/get_genotypes/get_subset_vcf.sh /srv/gsfs0/projects/montgomery/mohanas/africa/genotypes/21.minAC1.no_mask.with_related.id_filled.vcf.gz  /srv/gsfs0/projects/montgomery/mohanas/africa/genotypes/21.minAC1.no_mask.with_related.id_filled.599_samples  /srv/gsfs0/projects/montgomery/mohanas/africa/rna/sample_names.txt
```

> run_fastQTL.sh
```sh
#!/bin/bash


/srv/gsfs0/projects/montgomery/mohanas/tools/FastQTL-2.165.linux/bin/fastQTL.1.165.linux --vcf $1 --bed $2 --region $3 --out $4  --cov $5 --threshold 0.05 --include-samples samples.inc
```

```sh
qsub -cwd -M mohanas@stanford.edu -m ae -e EOs -o EOs -l h_vmem=200G /srv/gsfs0/projects/montgomery/mohanas/africa/rna/script/get_genotypes/add_unique_ids_fivePopulations.sh /srv/gsfs0/projects/montgomery/kameronr/Africa/genotype/VCFs/fivePopulations.21.minAC1.no_mask.with_related.recode.vcf  /srv/gsfs0/projects/montgomery/mohanas/africa/genotypes/fivePopulations.21.minAC1.no_mask.with_related.recode.id_filled.vcf /srv/gsfs0/projects/montgomery/mohanas/africa/genotypes/fivePopulations.21.minAC1.no_mask.with_related.recode.id_filled.vcf.gz
```
Per population fastqtl on durga
```
/srv/persistent/bliu2/tools/fastqtl/bin/fastQTL --vcf /srv/scratch/mohanas/africa/genotypes/21.minAC1.no_mask.with_related.id_filled.vcf.gz --bed /srv/scratch/mohanas/africa/rna/phenotypes/YRI_phenotypes.bed.gz --out EQTL-Results/21.minAC1.no_mask.with_related.id_filled.YRI_phenotypes.YRI_covariates_7 --cov /srv/scratch/mohanas/africa/rna/covariates/YRI_covariates_7.txt.gz --threshold 0.05 --region 21 --exclude-samples /srv/scratch/mohanas/africa/rna/exclude_samples/YRI.exc
```
Run fastQTL on durga when all the input files needed are ready:
>run_fastQTL.py
```py
import os
import glob
import time

PHENO_DIR_PATH = "/srv/scratch/mohanas/africa/rna/"
PHENO_FILE_SUFFIX = "phenotypes_mini.bed.gz" ## actual file = population | suffix

GENO_DIR_PATH = "/srv/scratch/mohanas/africa/genotypes/"
GENO_FILE_SUFFIX = ".minAC1.no_mask.with_related.id_filled.vcf.gz" ## actual file = chr # | suffix

COV_DIR_PATH = "/srv/scratch/mohanas/africa/rna/"
COV_FILES = ["covariates_34.txt.gz"]

RESULTS_DIR = "EQTL-Results/"

POPULATIONS = [""]
CHRS = ["21"]

NUM_CHUNKS_TOTAL = "500"

start = time.time()
ct = 0
for chrom in CHRS:
	for pop in POPULATIONS:
		for cov_file in COV_FILES:
			cmd1 = "/srv/persistent/bliu2/tools/fastqtl/bin/fastQTL --vcf " + GENO_DIR_PATH + chrom + GENO_FILE_SUFFIX + \
				   " --bed " + PHENO_DIR_PATH + pop + PHENO_FILE_SUFFIX + \
				   " --out " + RESULTS_DIR + chrom + GENO_FILE_SUFFIX[:-6] + pop + PHENO_FILE_SUFFIX[:-6] + cov_file[:-7] + \
				   " --cov " + COV_DIR_PATH + cov_file + \
				   " --threshold 0.05 --commands " + NUM_CHUNKS_TOTAL + " commands." + NUM_CHUNKS_TOTAL + ".txt" 

			print cmd1
			os.system(cmd1)

			fil = open("commands." + NUM_CHUNKS_TOTAL + ".txt")
			for line in fil:
				if line.strip().split()[-1].split(":")[0] == chrom:
					ct += 1
					print str(ct) + ") " + line.strip()
					os.system(line.strip())
					end = time.time()
					print "Time elapsed (in secs) = ", end-start



```

