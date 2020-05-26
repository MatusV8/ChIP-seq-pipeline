# ChIP-seq pipeline

This is a snakemake pipeline designed to analyse ChIP-seq projects in a clear and consistent way.

This pipeline can analyse both single-end (SE) and pair-end (PE) samples within a single ChIP-seq project.
The analysis includes:

** Merging multiple fastq files corresponding to a single sample
** Adapter trimming using and quality filtering using [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
** Quality control of sequencing using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
** Read mapping to a reference genome using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
** Marking duplicates using [samblaster](https://github.com/GregoryFaust/samblaster)
** Quality filtering of BAM files
** Quality control using [phantompeakqualtools](https://github.com/kundajelab/phantompeakqualtools) and [deeptools](https://deeptools.readthedocs.io/en/develop/)
** Peak calling using [MACS2](https://github.com/taoliu/MACS)
** [Irreproducibility Discovery Rate (IDR)](https://www.encodeproject.org/software/idr/) analysis of samples with at least two replicates
** Peak annotation using [HOMER](http://homer.ucsd.edu/homer/ngs/annotation.html)
** Creating a reporting summary using [multiQC](https://multiqc.info/)

## Prerequisites

To use this pipeline, install conda as described [here](https://bioconda.github.io/user/install.html).

Clone the pipeline into your current working directory

```
git clone https://github.com/MatusV8/ChIP-seq-pipeline.git
```
Create conda environment containing all required software for the pipeline and activate it.

```
conda env create -f ./env/ChIP-seq.yml
conda activate ChIP-seq
```

Check if the HOMER installation contains reference files for the genome you want to use using python3 `HOMER_config_check.py` script. If it genome is not installed, it will be automatically downloaded and installed.

```
# Example for mm10 genome
python3 ./bin/HOMER_config_check.py -g mm10
```

## Getting started

In order to analyse your set of ChIP-seq experiments, you need to define your experimental design and specify configuration for your analysis. This can be done by easily done by editing `Design.tsv` and `./config/config.yaml` files.

### Experimental design configuration

The experimental design is defined by editing Design.tsv file. Each row of corresponds to a single
sample and its control. Example of the Design.tsv is shown:

Cell_line | Condition | Factor | Replicate  | Strandeness | File | Control_file
-------- | --------  | -------- | --------  | -------- | -------- | --------
A        | 0h   | RNAPII | 1 | PE | RNAPII-A0h.fastq.gz | Control_A.fastq.gz
A        | 24h  | RNAPII | 1 | PE | RNAPII-A24h.fastq.gz | Control_A.fastq.gz
A        | 0h   | Oct4   | 1 | SE | Oct4-A0h.fastq.gz | Control_A.fastq.gz
A        | 24h  | Oct4   | 1 | SE | Oct4-A0h.fastq.gz  | Control_A.fastq.gz
B        | WT   | RNAPII | 1 | PE | RNAPII-WT1.fastq.gz  | Control_B.fastq.gz
B        | WT   | RNAPII | 2 | PE | RNAPII-WT2.fastq.gz | Control_B.fastq.gz
B        | WT   | Oct4   | 1 | SE | Oct4-WT1.fastq.gz | Control_B.fastq.gz
B        | WT   | Oct4   | 2 | SE | Oct4-WT2.fastq.gz | Control_B.fastq.gz


**Cell_line**, **condition**, **factor** and **replicate** columns specify experimental conditions of the ChIP sample, cell line used, experimental condition (treatment), immunoprecipitated factor and a replicate number. Replicate number have to be an integer (e.g. 1) and needs to be specified for each row. These three columns are used to create a unique sample ID. For example, `A-Oct4-24h-rep1` would uniquely identify a Oct4 ChIP-seq sample from the cell line A treated for 24 hours. Results are grouped by **factor** to make it easier to find.

**Strandeness** column accepts 'SE' or 'PE' values to determine the type of sequencing used.

**File** determines a path(s) of the sequencing file(s) 
in .fastq .sanfastq .fastq.gz or .sanfastq.gz formats. Multiple files can be added for samples sequenced 
in multiple lanes by separating paths by a comma. For example `./sample1-1_1.fastq.gz,./sample1-2_1.fastq.gz` would specify two files for the first read of the sample `sample1`. Please note that you don't have to specify path for the file containing second pairs. The software assume the presence of that file in the same directory as the first one. For example it looks for files `./sample1-1_2.fastq.gz` and `./sample1-2_2.fastq.gz`.

**Control_file** contains a path to the file used as a control for 
the sample. Multiple file paths can be specified as in **File**. A control file has to be specified to each ChIP sample. It is possible to use single control for multiple ChIP samples by simply specifying the same path for multiple ChIP samples. The pipeline will recognise this and prevents duplicated analysis of the same control.

### Configuration of the parameters

The `./config/config.yaml` config file contains various parameters used for the analysis. Before you start your analysis, you have to specify paths to **bowtie2_index** and **blacklist regions** (in a bed format), genome you want to use (e.g "mm10" or "hg18") and **effective genome size** for normalization. All parameters except **effective genome size** have to be quoted with `""`. Optionally, you can edit the default parameters for quality filtering, peak calling or figure generation. Please refer to the manuals listed [above]](#chip-seq-pipeline).

Pre-build **bowtie2 indexes** can be downloaded [here](http://bowtie-bio.sourceforge.net/manual.shtml).

**Effective genome sizes** can be found [here](https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html).


### Running the pipeline

Prior your analysis, you can check if there are no errors in the configuration using `-n` parameter (dry-run).
```
snakemake -j 10 -n
```
Note that `-j` has to be specified which tells snakemake how many threads are available to use.

When no errors are detected you can start the analysis by omitting the `-n` parameter.

```
snakemake -j 10
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
