# B-HIVE for single provirus transcriptomics

### 1. About
This repository contains all necessary files to map and compute HIV expression at single provirus resolution.

### 2. Prerequisites
The following packages are required to run the pipeline:
- Java Runtime environment 7
- Python 2.7
- Nextflow ([website](https://www.nextflow.io/))
- Docker ([website](https://docker.com/))

#### Installing prerequisites:
Docker and Nextflow must be downloaded following the instructions on their official websites.
To install *JRE7* and *Python 2.7* in *Ubuntu 14.04*, type in a terminal:
```
sudo apt-get install openjdk-7-jre python2.7
```

### 3. Running the pipeline:
#### 3.1. Download sources
Download and extract all files from this repository. Use either `git`:
```
git clone https://github.com/gui11aume/BHIVE_for_single_provirus_transcriptomics/
```
or `wget` if you don't have `git` installed in your computer:
```
wget https://github.com/gui11aume/BHIVE_for_single_provirus_transcriptomics/archive/master.zip
unzip master.zip
rm master.zip
mv BHIVE_for_single_provirus_transcriptomics-main BHIVE_for_single_provirus_transcriptomics
```
then `cd` to the source directory:
```
cd BHIVE_for_single_provirus_transcriptomics
```
#### 3.2. Edit configuration files
You must edit `map.cfg`. `expr.cfg` and (optional) `demux.cfg`, to let the pipeline know the path to the
sequencing files that will be processed. Detailed instructions are provided inside the files.

In this example, all `.cfg` files are set to work with the example files provided along with the pipeline
sources (see `examples` folder).

#### 3.3. Generate a bwa index (optional)
Run this step only if you don't have yet an index of the Human Genome reference for `bwa`
([repository](https://github.com/lh3/bwa)). To build a `bwa` index of *hg19*:
```
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/chromFa.tar.gz
tar -Oxf chromFa.tar.gz > ref-genome.fasta
rm chromFa.tar.gz
docker run -v $(pwd):/data ezorita/bioinformatics bwa index /data/ref-genome.fasta
```

#### 3.4. Map HIV integrations
Once the paths to the mapping-iPCR files have been specified in `map.cfg`, run the mapping pipeline:
```
nextflow run map.nf --index ref-genome.fasta -with-docker ezorita/bioinformatics
```

where `ref-genome.fasta` is the path to the bwa-index. After a successful run, it will 
generate a file named `hiv_integrations.txt` with the following columns:
- **brcd:**   provirus barcode.
- **chr:**    integration chromosome.
- **locus:**  integration nucleotide.
- **strand:** integration strand.
- **reads:**  count of provirus molecules sequenced.
- **mapq:**   confidence score (from 0 to 100).
- **rep:**    replicate number.

#### 3.5. Demultiplex barcode-PCR files (optional)
If the DNA and RNA barcodes were multiplexed in the same sequencing pool using DNA indices, an additional step
is required to demultiplex the barcodes and assign them to the correct experiment replicate. The files and
indices should be indicated in `demux.cfg` (detailed instructions inside the file). Then run:
```
python demux.py demux.cfg
```
The new demultiplexed files will be generated in FASTQ format, one for each DNA/RNA replicate.

#### 3.6. Compute provirus expression
Finally, the paths to the demultiplexed DNA/RNA files should be specified in `expr.cfg`. To run the expression pipeline, type:
```
nextflow run expr.nf --integs hiv_integrations.txt -with-docker ezorita/bioinformatics
```
where `hiv_integrations.txt` is the path to the file generated in the mapping step. The output of this process is a file
named `hiv_expression.txt` with the following columns:
- **brcd:**   provirus barcode.
- **rep:**    replicate number.
- **chr:**    integration chromosome.
- **locus:**  integration nucleotide.
- **strand:** integration strand.
- **reads:**  count of provirus molecules sequenced.
- **mapq:**   confidence score (from 0 to 100).
- **dna:**    count of DNA barcode molecules sequenced.
- **rna:**    count of mRNA barcode molecules sequenced.
- **exprscore:** balanced expression score.

Also, a scatter plot figure is generated for each technical replicate (see `figures/` folder).

### 4. References


