# DupCaller
[![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://osf.io/s93d5/wiki/home/) [![License](https://img.shields.io/badge/License-BSD\%202--Clause-orange.svg)](https://opensource.org/licenses/BSD-2-Clause) [![Build Status](https://travis-ci.com/AlexandrovLab/SigProfilerMatrixGenerator.svg?branch=master)](https://app.travis-ci.com/AlexandrovLab/SigProfilerMatrixGenerator)
[![Uptime Robot status](https://img.shields.io/uptimerobot/status/m795312784-02766a79f207f67626cef289)](https://stats.uptimerobot.com/jjqW4Ulymx)

# DupCaller
DupCaller is a universal tool for calling somatic mutations and calculating somatic mutational burden from barcoded error-corrected next generation sequencing (ecNGS) data with matched normal (e.x. NanoSeq). 

**Prerequisites***
DupCaller requires python>=3.10 to run. Earlier versions may be sufficient to run DupCaller and have not been tested. 
The complete DupCaller pipeline also requires the following tools for data preprocessing. The versions are used by the developer and other versions may or may not work.
  * BWA          version 0.7.17 (https://bio-bwa.sourceforge.net)
  * GATK         version 4.2.6 (https://github.com/broadinstitute/gatk/releases)
  * mosdepth     version 0.3.3 (https://github.com/brentp/mosdepth)
**INSTALLATION**
The tool uses pip for installing scripts and prerequisites. To install DupCaller, simply clone this repository and install via pip:
    '''bash
    git clone https://github.com/AlexandrovLab/DupCaller.git
    cd DupCaller
    pip install .
    '''

**Pipeline**

#### Trim barcodes from reads:
DupCallerTrim.py is a scripts that can extract 5-prime barcodes from paired-end fastqs. The usage is as follows:
    '''bash
    DupCallerTrim.py -i read1.fq -i2 read2.fq -p barcode_pattern -o sample_name
    '''
where
'read1.fq' and 'read2.fq' are fastq files from read1 and read2 of the paired-end sequencing data, respectively. Both unzipped and gzip compressed files can be correctly processed.
'barcode_pattern' is the pattern of barcodes starting from the 5-prime end, with N representing a barcode base and X representing a skipped base. The notation is similar to what has been used in UMI-tools(https://github.com/CGATOxford/UMI-tools). For example, NanoSeq uses 3-base barcode followed by 4 constant bases, so the pattern should be NNNXXXX.
'sample_name' is the prefix of output paired fastqs. After run complete, '{sample_name}_1.fastq' and '{sample_name}_2.fastq' will be generated.
The barcodes will be recorded in the each read name as {original_read_name}:{read1_barcode}+{read2_barcode}

If the matched normal is prepared in the same way as the sample, also apply the trimming with the same scheme to the matched normal fastqs. For traditional bulk normal, trimming is not needed.

#### Align trimmed fastqs
Use a DNA NGS aligner, such as BWA-MEM, to align the trimmed fastqs of both sample and matched normal from the last step. Notice that GATK requires read group ID,SM and PL to be set, so adding those tags during bwa alignment is recommended. For example:
    '''bash
    bwa mem -t {threads} -R "@RG\tID:{sample_name}\tSM:{sample_name}\tPL:ILLUMINA" reference.fa {sample_name}_1.fastq {sample_name}_2.fastq | samtools sort -@ {threads} > {sample_name}.bam
    samtools index -@ {threads} {sample_name}.bam
    '''
where
'threads' is the number of cores used for aligning
'reference.fa' is the reference genome fasta file
'{sample_name}_1.fastq' and '{sample_name}_2.fastq' are trimmed fastq file from last step.

#### (Optional) Running BQSR on sample and matched normal
Base quality score recalibration (BQSR) of the bam file is optional and can be performed following https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-. The effect of BQSR on ecNGS variant calling is theoretically beneficial, since accurate base quality score can improve variant calling accuracy of lowly-duplicated read groups. 
 1. generate BQSR reports with GATK BaseRecalibrator:https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator
 2. Apply BQSR on the bam with GATK applyBQSR: https://gatk.broadinstitute.org/hc/en-us/articles/360036856671-ApplyBQSR

#### MarkDuplicates with optical duplicates tags and new read name configuration
Run GATK MarkDuplicates on sample and matched-normal bams. Notice that optical duplicates and PCR duplicates should be treated differently in ecNGS variant calling, so the TAGGING_POLICY of GATK MarkDuplicates should be set to OpticalOnly to differentiate optical duplicate from PCR duplicate. Also, since the read name of trimmed fastq is non-traditional, the READ_NAME_REGEX option should also be set to "(?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*$". The MarkDuplicates commands should be looking like this:
'''bash
gatk MarkDuplicates -I sample.bam -O sample.mkdped.bam -M sample.mkdp_metrics.txt --READ_NAME_REGEX "(?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*$" --TAGGING_POLICY OpticalOnly
'''

#### Running mosdepth on normal bam file
We take advantage of the per-base coverage file generated by mosdepth 

**PARAMETERS**
| Category | Parameter | Variable Type | Parameter Description |
| ------ | ----------- | ----------- | ----------- |
| Required |  |  |  |
|  | project | String | The name of the project. |
|  | reference_genome | String | The name of the reference genome. Full list of genomes under **Supported Genomes** section. Supported values include the following: {c_elegans, dog, ebv, GRCh37, GRCh38, mm9, mm10, mm39, rn6, yeast} |
|  | path_to_input_files | String | The path to the input files. |
| Optional |  |  |  |
|  | exome | Boolean | Downsamples mutational matrices to the exome regions of the genome. Default value False. |
|  | bed_file | String | Downsamples mutational matrices to custom regions of the genome. Requires the full path to the BED file. Default value None. Note: BED file header is required. |
|  | chrom_based | Boolean | Outputs chromosome-based matrices. Default value False. |
|  | plot | Boolean | Integrates with SigProfilerPlotting to output all available visualizations for each matrix. Default value False. |
|  | tsb_stat | Boolean | Outputs the results of a transcriptional strand bias test for the respective matrices. Default value False. |
|  | seqInfo | Boolean | Ouputs original mutations into a text file that contains the SigProfilerMatrixGenerator classificaiton for each mutation. Default value True. |
|  | cushion | Integer | Adds an Xbp cushion to the exome/bed_file ranges for downsampling the mutations. Default value 100. |
|  | volume | String | Path to SigProfilerMatrixGenerator's volume where reference genomes will be saved and loaded from. Useful for docker installations. Default is None. |


**INPUT FILE FORMAT**

This tool currently supports maf, vcf, simple text file, and ICGC formats. The user must provide variant data adhering to one of these four formats. If the user’s files are in vcf format, each sample must be saved as a separate files.


**OUTPUT FILE STRUCTURE**

The output structure is divided into three folders: input, output, and logs. The input folder contains copies of the user-provided input files. The outputfolder contains
a DBS, SBS, ID, and TSB folder (there will also be a plots folder if this parameter is chosen). The matrices are saved into the appropriate folders. The logs folder contains the error and log files for the submitted job.

## STRUCTURAL VARIANT MATRIX GENERATION

### INPUT FORMAT:

***First six columns are required, and either the column "svclass" (deletion, translocation, tandem-duplication, or inversion) or the columns "strand1" & "strand2" (BRASS convention) must also be present***


### Example with SV class present (tsv or csv file):

| chrom1 | start1 | end1 | chrom2 | start2 | end2 | svclass |
| :-----: | :-: | :-: | :-: | :-: | :-: | :-: |
| 19 | 21268384 | 21268385 | 19 | 21327858 | 21327859 | deletion

### Example without SV class present (tsv or csv file):

| chrom1 | start1 | end1 | chrom2 | start2 | end2 | strand1 | strand2
| :-----: | :-: | :-: | :-: | :-: | :-: | :-: | :-: |
| 19 | 21268384 | 21268385 | 19 | 21327858 | 21327859 | + | +

### Quick Start Example: ###

```
#navigate to SVMatrixGenerator directory and start python3 interpreter

from SigProfilerMatrixGenerator.scripts import SVMatrixGenerator as sv
input_dir = "./SigProfilerMatrixGenerator/references/SV/example_input/560-Breast" #directory which contains collection of bedpe files (one per sample)
output_dir = "./SigProfilerMatrixGenerator/references/SV/"
project = "560-Breast"
sv.generateSVMatrix(input_dir, project, output_dir)
```
**Alternatively, you can run directly from the command line:**
```
python3 ./SigProfilerMatrixGenerator/scripts/SVMatrixGenerator.py ./SigProfilerMatrixGenerator/references/SV/example_input/560-Breast 560-Breast ./SigProfilerMatrixGenerator/references/SV/example_output/ #provide input_dir, project, output_dir as command-line arguments
```
## OUTPUT:
1. Annotated bedpe file - a file with each SV annotated with its type, size bin, and clustered/non-clustered status
2. Aggregate SV plot - a summary plot showing the average number of events in each channel for the whole cohort of samples
3. SV Matrix - a 32 X n matrix (where n is the number of samples) that can be used to perform signature decomposition, clustering, etc.

## COPY NUMBER MATRIX GENERATION

In order to generate a copy number matrix, provide the an absolute path to a multi-sample segmentation file obtained from one of the following copy number calling tools (if you have individual sample files, please combine them into one file with the first column corresponding to the sample name):

1. ASCAT
2. ASCAT_NGS
3. SEQUENZA
4. ABSOLUTE
5. BATTENBERG
6. FACETS
7. PURPLE
8. TCGA

In addition, provide the name of the project and the output directory for the resulting matrix. The final matrix will be placed in a folder with the name of the project in the directory specified by the output path.

**An example to generate the CNV matrix is as follows:**

$ python3
```
>>from SigProfilerMatrixGenerator.scripts import CNVMatrixGenerator as scna
>>file_type = "BATTENBERG"
>>input_file = "./SigProfilerMatrixGenerator/references/CNV/example_input/Battenberg_test.tsv" #example input file for testing
>>output_path = "/Users/azhark/iCloud/dev/CNVMatrixGenerator/example_output/"
>>project = "Battenberg_test"
>>scna.generateCNVMatrix(file_type, input_file, project, output_path)

```

**Alternatively, you can run directly from the command line:**

```
python ./SigProfilerMatrixGenerator/scripts/CNVMatrixGenerator.py BATTENBERG ./SigProfilerMatrixGenerator/references/CNV/example_input/Battenberg_test.tsv BATTENBERG-TEST ./SigProfilerMatrixGenerator/references/CNV/example_output/
```

**SUPPORTED GENOMES**

This tool currently supports the following genomes:

GRCh38.p12 [GRCh38] (Genome Reference Consortium Human Reference 38), INSDC
Assembly GCA_000001405.27, Dec 2013. Released July 2014. Last updated January 2018. This genome was downloaded from ENSEMBL database version 93.38.

GRCh37.p13 [GRCh37] (Genome Reference Consortium Human Reference 37), INSDC
Assembly GCA_000001405.14, Feb 2009. Released April 2011. Last updated September 2013. This genome was downloaded from ENSEMBL database version 93.37.

GRCm39 [mm39] (Genome Reference Consortium Mouse Reference 39), INSDC
Assembly GCA_000001635.9, Jun 2020. Last updated August 2020. This genome was downloaded from ENSEMBL database version 103.

GRCm38.p6 [mm10] (Genome Reference Consortium Mouse Reference 38), INDSDC
Assembly GCA_000001635.8, Jan 2012. Released July 2012. Last updated March 2018. This genome was downloaded from ENSEMBL database version 93.38.

GRCm37 [mm9] (Release 67, NCBIM37), INDSDC Assembly GCA_000001635.18.
Released Jan 2011. Last updated March 2012. This genome was downloaded from ENSEMBL database version release 67.

Rnor_6.0 [rn6] INSDC Assembly GCA_000001895.4, Jul 2014. Released Jun 2015. Last updated Jan 2017.
This genome was downloaded from ENSEMBL database version 96.6.

Epstein-Barr Virus [EBV] NC_007605.1, Nov 2005. Last updated Aug 2018. This genome was downloaded from the NCBI database: https://www.ncbi.nlm.nih.gov/nuccore/82503188/.

CanFam3.1 [dog] GCA_000002285.2, Sep 2011. Last updated Jun 2019. This genome was downloaded from ENSEMBL database version 100.

WBcel235 [c_elegans] GCA_000002985.3, Oct 2014. Last updated Jan 2019. This genome was downloaded from ENSEMBL database version 100.

*One can specify "_havana" to the end of the genome to include annotations in t-cell receptor genes and IG clusters (available for GRCh37, GRCh38, and mm10).

**LOG FILES**

All errors and progress checkpoints are saved into *sigProfilerMatrixGenerator_[project]_[genome].err* and *sigProfilerMatrixGenerator_[project]_[genome].out*, respectively.
For all errors, please email the error and progress log files to the primary contact under CONTACT INFORMATION.

**UNIT TESTS**

Unit tests can be run with the following commands:

```bash
python setup.py sdist
pip install .[tests]
pytest tests
```

**END-TO-END TESTS**

An integration test can be run with the following commands:

```bash
wget ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerMatrixGenerator/GRCh37.tar.gz -P ./src/
pip install .
SigProfilerMatrixGenerator install GRCh37
python3 test.py -t GRCh37
```

**CITATION**

For SBSs, DBSs, and INDELs, please cite the following paper:

Bergstrom EN, Huang MN, Mahto U, Barnes M, Stratton MR, Rozen SG, and Alexandrov LB (2019) SigProfilerMatrixGenerator: a tool for visualizing and exploring patterns of small mutational events. **BMC Genomics** 20, Article number: 685.
https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-6041-2

For SVs and CNVs, please cite the following paper:

Khandekar A, Vangara R, Barnes M, Díaz-Gay M, Abbasi A, Bergstrom EN, Steele CD, Pillay N, and Alexandrov LB (2023) Visualizing and exploring patterns of large mutational events with SigProfilerMatrixGenerator. **BMC Genomics** 24, Article number: 469.
https://doi.org/10.1186/s12864-023-09584-y

**COPYRIGHT**

Copyright (c) 2019, Erik Bergstrom [Alexandrov Lab] All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

**CONTACT INFORMATION**

Please address any queries or bug reports to Erik Bergstrom at ebergstr@eng.ucsd.edu. Please address any queries or bug reports related to CNV's or SV's to Azhar Khandekar at akhandek@eng.ucsd.edu. Additional support can be provided by Mark Barnes at mdbarnes@health.ucsd.edu.