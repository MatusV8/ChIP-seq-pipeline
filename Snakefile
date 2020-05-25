##########################################################
#
#------------ ChIP-seq Snakemake Pipeline ---------------#
#
##########################################################
# Author: Matus Vojtek

#-----------------------------------------
# Setup
#-----------------------------------------

import pandas as pd
import os
import re
import sys
import re
from pathlib import Path


# List of standard chromosomes excluding chrM in different species. Keys are UCSC assemblies
chromosomes = {
    "mm9": ["chr" + str(i) for i in range(1,20)] + ["chrX", "chrY"],
    "mm10": ["chr" + str(i) for i in range(1,20)] + ["chrX", "chrY"],
    "hg38": ["chr" + str(i) for i in range(1,23)] + ["chrX", "chrY"],
    "hg19": ["chr" + str(i) for i in range(1,23)] + ["chrX", "chrY"],
    "dm6": ["chr4", "chr2", "chr2L", "chr2R", "chr3L", "chr3M", "chrY", "chrX"],
    "dm3": ["chr4", "chr2", "chr2L", "chr2R", "chr3L", "chr3M", "chrY", "chrX"],
    "danre10": ["chr" + str(i) for i in range(1,26)] + ["chrX", "chrY"],
    "danre11": ["chr" + str(i) for i in range(1,26)] + ["chrX", "chrY"] 
}


configfile: "config/config.yaml"

# Load in the design table
design = pd.read_csv('Design.tsv', sep = '\t')

#-------------------------------------------------
# Pre-process design table and check requirements
#-------------------------------------------------

# Check the Strandeness is in the correct format
if set(design.Strandeness) != {'SE', 'PE'}:
    sys.exit("Strandeness in the design.csv has a wrong format! Please use 'SE' or 'PE'.")

# Split File column by comma to get list of fastq files per sample
design.File = design.File.str.split(',')
design.Control_file = design.Control_file.str.split(',')

# Create name for each sample based on design table
# Multiple replicates can have the same name
design['Name'] = design.Cell_line + '-' + design.Factor + '-' + design.Condition

# Check if the ID names exist
if 'ID' in design.columns:
    
    print('ID names found. IDs used for  naming.')
    
    # Sort design dataframe by column
    design = design.sort_values(by=['ID'])
    
else:
    
    print('ID names not found.\n\nCreating ID names....\n')
    
    # Create unique ID for each sample which is human readable.
    # ID looks like Cell_line-Transcription_Factor-Condition-rep_N 
    design['ID'] = design.Cell_line + '-' + design.Factor + '-' + design.Condition + '-rep' + design.Replicate.apply(str)
    
    # Sort design dataframe by colum which is important later when dealing with inputs
    design = design.sort_values(by=['ID'])
    
# Check if the ID names are unique. If not, shut down with error message
if len(design.ID) != len(set(design.ID)):
    sys.exit("ID names are not unique !!!\nChange ID names to be unique")

# Check if Control_ID exists
if 'Control_ID' in design.columns:
    print('Control_ID names found. Control_ID used for naming.')
    
    # Add sample name column for controls
    design['Control_name'] = design['Name']
    
    # Because more than 1 sample can use the same input.
    # Check which controls are duplicated and replace them with the first one
    for ind in design.index:
        indexes = [i == design.Control_file[ind] for i in design.Control_file]
        design.loc[indexes, 'Control_ID'] = design.loc[indexes, 'Control_ID'].tolist()[0]
        design.loc[indexes, 'Control_name'] = design.loc[indexes, 'Control_name'].tolist()[0]

else:
    print('Creating Control_IDs ....\n')
    design['Control_ID'] = design.Name + '-Input-rep' + design.Replicate.astype('str')
    design['Control_name'] = design['Name']
    # Because more than 1 sample can use the same input.
    # Check which are duplicated and replace them with the first one
    for ind in design.index:
        indexes = [i == design.Control_file[ind] for i in design.Control_file]
        design.loc[indexes, 'Control_ID'] = design.loc[indexes, 'Control_ID'].tolist()[0]
        design.loc[indexes, 'Control_name'] = design.loc[indexes, 'Control_name'].tolist()[0]

# Create Files Dataframe containing both Control and IP files
# This dataframe is used to create targets for individual fastq files
Files = pd.concat([design.loc[:,['ID', 'File' ,'Strandeness', 'Name']],
        design.loc[:,['Control_ID', 'Control_file' ,'Strandeness', 'Name']].rename(columns = {'Control_ID' : 'ID', 'Control_file' : 'File'})],
    ignore_index = True)

# Remove duplicated files if exist
# might be redundant
Files = Files.drop_duplicates(subset = 'ID')

# Check the paths of the files
File_paths = []
[File_paths + i for i in Files.File]

for path in File_paths:
    if not os.path.exists(path):
        sys.exit(path + " not found !!!")
    

#---------------------------
# Determine file extensions
#---------------------------

# Add file extensions
extensions = ['.fq', '.fastq', '.sanfastq', '.fq.gz', '.fastq.gz', 'sanfastq.gz'] 

# Loop over File dataframe to check if the
# files are in the correct format and append it
# as a "Suffix" column to the Files dataframe
for i in Files.index:
    for ext in extensions:
       ext_found = [x.endswith(ext) for x in Files.File[i]]
       
       # If extension found in all files then report it
       if all(ext_found):
            extension = ext
            
       # If not all files have the same extension report error
       elif any(ext_found):
           ErrorMsg = 'Different file formats supplied for sample ' + Files.ID[i]
           sys.exit(ErrorMsg) 
   
    # if Extension found add it to the table
    if extension:
        Files.loc[i ,'Suffix'] = extension
    else:
        ErrorMsg = 'File format is not recognised. Options are .fq .fq.gz .fastq or .fastq.gz'
        sys.exit(ErrorMsg)
        
        
#------------------------------------------------------------
# Add second read files for PE samples to the File dataframe
#------------------------------------------------------------

# Create new data frame with all files. This expands the original data frame for the second strand files
Files1 = pd.DataFrame(columns = ['ID', 'File',])
for i in Files.index:
    nfiles = len(Files.loc[i, 'File'])
    suffix = Files.loc[i, 'Suffix']
    path1 = Files.loc[i, 'File']
    
    # If a sample is 'PE', automatically create second read sample
    # with the same path
    if Files.loc[i, 'Strandeness'] == 'PE':
       names =  [Files.loc[i, 'ID'] + str(j) for j in ['_1', '_2']]
       suff1 = '_1' + suffix
       suff2 = '_2' + suffix
       
       # Because there might be multiple fastq files for single sample
       # stored in a list, loop over the list to create second read file paths
       path2 = [x.replace(suff1, suff2) for x in path1]

       Files1 = Files1.append({'ID': names[0], 'File': path1}, ignore_index=True)
       Files1 = Files1.append({'ID': names[1], 'File': path2}, ignore_index=True)
    
    
    else:
        name = Files.loc[i, 'ID']
        Files1 = Files1.append({'ID': name, 'File': path1}, ignore_index=True)

SAMPLES = pd.DataFrame({'Name': design.Name.drop_duplicates()})
SAMPLES['N_replicates'] = [design.ID.str.contains(x).sum() for x in SAMPLES.Name]

#-------------------------------
# Create targets for BAM files 
#------------------------------

# Create data frame for all bam files including pseudoreplicates and merged files
IDR_targets = pd.DataFrame( columns = ['Name', 'ID', 'Control_ID', 'Control_name'])

for name in design.Name.drop_duplicates():
    replicates = design.loc[design.Name == name, 'ID']                  # Get sample IDs 
    controls = design.loc[design.Name == name, 'Control_ID']
    names = design.loc[design.Name == name, 'Name']
    Controlnames = design.loc[design.Name == name, 'Control_name']
    controlname = Controlnames.drop_duplicates()
    nreplicates = len(replicates)                 # Get number of replicates for the sampple
    
    if nreplicates > 1:
        PR1 = replicates + '_PR1'
        PR2 = replicates + '_PR2'
        PR = replicates.append([PR1,PR2])
        
        PR1input = controls + '_PR1'
        PR2input = controls + '_PR2'
        PRinput = controls.append([PR1input,PR2input])
        
        # add pseudoreplicates to the table
        PR_pd = pd.concat([pd.concat([names] * 3), 
                         PR, PRinput,
                         pd.concat([Controlnames] * 3)], axis = 1).reset_index(drop =True)
        # Merged file                                                 # If the sample is replicate,
        mergedName = name + '-merged'                                 # add it to the target list,
        mergedPR1 = mergedName + '_PR1'
        mergedPR2 = mergedName + '_PR2'
        mergedPR = pd.Series([mergedName, mergedPR1 ,mergedPR2], name = 'ID')
        
        if len(controls.drop_duplicates()) == len(replicates):
            
            mergedInput = controlname + '-Input-merged'
            mergedInputPR1 = mergedInput + '_PR1'
            mergedInputPR2 = mergedInput + '_PR2'
            mergedInputPR = pd.concat([mergedInput, mergedInputPR1 ,mergedInputPR2]).reset_index(drop =True)
        
        else:
            print("Number of replicates for sample ", name, " and its control is different!")
            mergedInput = controls.drop_duplicates()
            mergedInputPR1 = mergedInput + '_PR1'
            mergedInputPR2 = mergedInput + '_PR2'
            mergedInputPR = pd.concat([mergedInput, mergedInputPR1 ,mergedInputPR2]).reset_index(drop =True)
        
        # add merged files to the table
        mergedPR_pd = pd.concat([ pd.Series([name]* 3),mergedPR, mergedInputPR,
                         pd.concat([controlname] * 3).reset_index(drop =True)],
                         axis = 1, ignore_index = True)
        mergedPR_pd.set_axis(['Name', 'ID', 'Control_ID', 'Control_name'], axis = 1, inplace=True)
        # append the tables to the final table
        IDR_targets = pd.concat([IDR_targets, PR_pd, mergedPR_pd], ignore_index = True)

BAMTARGETS = IDR_targets.append(design.loc[:,['Name', 'ID', 'Control_ID', 'Control_name']]).drop_duplicates().reset_index(drop = True)

#-----------------------
# Create multiQC header 
#-----------------------

multiqc_header = {'genome': config['genome'],
                  'Number of samples': len(SAMPLES.Name)}

#-----------------------------------
# Target rule to create all targets
#-----------------------------------

rule all:
    input:
        list('01_Fastq/Raw/' + Files1.ID + '.fastq.gz'),
        list('01_Fastq/Trimmed/' + Files1.ID + '.fastq.gz'),
        list('01_Fastq/FastQC/' + Files1.ID + '_fastqc.zip'),
        list('01_Fastq/FastQC/' + Files1.ID + '_fastqc.html'),
        list('00_Logs/Flagstat/Before/' + Files.ID + '.txt'),
        list('00_Logs/Flagstat/After/' + Files.ID + '.txt'),
        list('02_Bam/' + BAMTARGETS.Name + '/' + BAMTARGETS.ID + '.bam'),
        list('02_Bam/' + BAMTARGETS.Name[~BAMTARGETS.ID.str.contains('_PR')] + 
        '/' + BAMTARGETS.ID[~BAMTARGETS.ID.str.contains('_PR')] + '.bam.bai'),
        list('03_BigWigs/Coverage/' + Files.ID + '.bw'),
        list('03_BigWigs/BigWigSummary/' + design.Factor.drop_duplicates() + '.npz'),
        list('03_BigWigs/Input_normalised/' + design.ID + '.InputNormalised.bw'),
        list('04_QC/FingerPrint/' + design.Factor.drop_duplicates() + '.pdf'),
        list('04_QC/FingerPrint/' + design.Factor.drop_duplicates() + '.qualmat.tab'),
        list('04_QC/FingerPrint/' + design.Factor.drop_duplicates() + '.counts.tab'),
        list('04_QC/PCA/' + design.Factor.drop_duplicates() + '.pdf'),
        list('04_QC/PCA/' + design.Factor.drop_duplicates() + '.dat'),
        "04_QC/PCA/all.dat",
        list('04_QC/Cross_Correlation/' + Files.Name + '/' + Files.ID + '.cc.qc'),
        list('04_QC/Correlation/' + design.Factor.drop_duplicates() + '.pdf'),
        list('04_QC/Correlation/' + design.Factor.drop_duplicates() + '.tsv'),
        "04_QC/Correlation/all.tsv",
        list('05_Peaks/' + design.Name + '/' + design.ID + '_peaks.narrowPeak'),
        list('05_Peaks/' + IDR_targets.Name + '/IDR/' + IDR_targets.ID + '_peaks.narrowPeak'),
        "04_QC/QC_summary.csv", # FRIP analysis
        "04_QC/IDR_summary.tsv", # IDR analysis
        "10_multiQC/multiQC_log.html"

#----------------------------------------------
# Merge multiple fastq files for single sample
#----------------------------------------------

rule merge_fastq:
    input: lambda wildcards: (Files1.File[Files1.ID == wildcards.sample].tolist()[0])
    output: "01_Fastq/Raw/{sample}.fastq.gz"
    run: 
        if len(input) > 1:
            shell(
                r"""cat {input} > {output}"""
                )
        else:
            shell(r"""
                ln {input} {output}
                """)

#----------------------------------
# Trim and quality filter PE reads
#----------------------------------

rule trim_galore_pe:
    input: "01_Fastq/Raw/{sample}_1.fastq.gz", '01_Fastq/Raw/{sample}_2.fastq.gz'
    output:
        out1 = "01_Fastq/Trimmed/{sample}_1.fastq.gz",
        out2 = "01_Fastq/Trimmed/{sample}_2.fastq.gz",
    threads: 4
    log:
        "00_Logs/trim_galore/{sample}.log"
        
    shell:
        r"""
        trim_galore -q 30 -o 01_Fastq/Trimmed/ \
        --paired \
        --basename {wildcards.sample} \
        -j {threads} \
        {input} 2> {log}
        
        # Rename files
        mv 01_Fastq/Trimmed/{wildcards.sample}_val_1.fq.gz {output.out1}
        mv 01_Fastq/Trimmed/{wildcards.sample}_val_2.fq.gz {output.out2}
        """

#----------------------------------
# Trim and quality filter SE reads
#----------------------------------

rule trim_galore_se:
    input: "01_Fastq/Raw/{sample}.fastq.gz"
    output:
        out = "01_Fastq/Trimmed/{sample}.fastq.gz",
        report = "01_Fastq/Trimmed/{sample}.fastq.gz_trimming_report.txt"
    wildcard_constraints: sample = '.+(?<!_\d)'
    threads: 4
    log:
         "00_Logs/trim_galore/{sample}.log"
    shell:
        r"""
        trim_galore -q {config[trim_galore][q]} -o 01_Fastq/Trimmed/ \
        -j {threads} \
        --basename {wildcards.sample} \
        {input} 2> {log}
        
        # Rename files
        mv 01_Fastq/Trimmed/{wildcards.sample}_trimmed.fq.gz {output.out}
        """
        
#----------------------------------
# QC of fastq files after trimming
#----------------------------------

rule fastqc:
    input:
        "01_Fastq/Trimmed/{sample}.fastq.gz"
    output:
        html="01_Fastq/FastQC/{sample}_fastqc.html",
        zip="01_Fastq/FastQC/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    
    log:
        "00_Logs/FastQC/{sample}.log"
    shell:
        r"""
        fastqc {input} -o 01_Fastq/FastQC/ 2> {log}
        """

#-----------------------------------------
# Genome alignment and duplicates marking
#-----------------------------------------

# Function to determine input files for bowtie2 based on
# sequencing (PE / SE)

def align_strand(wc):
    nme = wc.sample + wc.ext
    strand = Files.loc[Files['ID'] == nme, 'Strandeness'].values[0]
    if strand == "PE":
        return(['01_Fastq/Trimmed/{sample}{ext}_1.fastq.gz',
        '01_Fastq/Trimmed/{sample}{ext}_2.fastq.gz'])
    else:
        return('01_Fastq/Trimmed/{sample}{ext}.fastq.gz')

rule align:
    output: 
        bam = temp("02_Bam/{sample}/{sample}{ext}-Raw.bam")
    input: 
        align_strand
    log:
        bw = "00_Logs/Bowtie2/{sample}{ext}.txt",
        samblaster = "00_Logs/Samblaster/{sample}{ext}.txt"
    wildcard_constraints: ext = '-Input-rep\d|-rep\d'
    threads: 20
    run:
        # If there are two fastq for single file => PE
        if len(input) == 2:
            shell(
            r"""
            bowtie2 -x {config[bw2_index]} -p {threads} \
            -1 {input[0]} \
            -2 {input[1]} 2> {log.bw} | samblaster 2> {log.samblaster} | samtools sort > {output.bam}
            """)
            
        # Else it is a SE sample
        else:
            shell(
            r"""
            bowtie2 -x {config[bw2_index]} -p {threads} {input} \
            2> {log.bw} | samblaster --ignoreUnmated 2> {log.samblaster} | samtools sort > {output.bam}
            """)

#-----------------------------------------
# Flagstat QC before read filtering
#-----------------------------------------

rule flagstat_before:
    output: "00_Logs/Flagstat/Before/{sample}{ext}.txt"
    input: rules.align.output.bam
    wildcard_constraints: ext = '-Input-rep\d|-rep\d'
    shell:
        """samtools flagstat {input} > {output} """

#----------------------------------------------
# Remove reads with improper alignment, pairing
# or marked as duplicate
#----------------------------------------------

rule filter_bam:
    output:
        # BAM files containing all reads will be deleted
        # after filtered BAM file is created
        tmpout = temp("02_Bam/{sample}/{sample}{ext}.tmp.bam"),
        tmpindex = temp("02_Bam/{sample}/{sample}{ext}.tmp.bam.bai"),
        out = protected("02_Bam/{sample}/{sample}{ext}.bam")
    input:  rules.align.output.bam
    wildcard_constraints: ext = '-Input-rep\d|-rep\d'
    params: 
        flag = lambda wc: 1804 if (Files.loc[Files.ID == wc.sample + wc.ext]['Strandeness'].values[0] == 'PE') else 1796,
        chr = chromosomes[ config['genome'] ]
    shell:
        r"""
        samtools view -b {input} -F {params.flag} \
        | bedtools intersect -v -abam stdin -b {config[blacklist]} > {output.tmpout}
        
        # Index bam file
        samtools index {output.tmpout} {output.tmpindex}
        
        # Remove reads mapped to non-standard chromosomes and chrM
        samtools view -b {output.tmpout} {params.chr} > {output.out}
        """

# Function to determine folder for a sample based on Name

def group_lookup(wc):
    grp = Files.loc[Files.ID == wc.sample, 'Name']
    return("02_Bam/" + grp + "/{sample}.bam")
    
#-----------------------------------------
# Flagstat QC after read filtering
#-----------------------------------------

rule flagstat_after:
    output: "00_Logs/Flagstat/After/{sample}.txt"
    input: group_lookup
    wildcard_constraints: sample = '|'.join(SAMPLES.Name.tolist())
    shell:
        """samtools flagstat {input} > {output} """

# Function to create inputs to merge replicates
def n_replicates(wc):
    nreplicates = SAMPLES.loc[SAMPLES.Name == wc.grp, 'N_replicates'].values[0]
    return(['02_Bam/{grp}/{grp}-rep' + str(i) + '.bam' for i in range(1, nreplicates+ 1)])
    
#-----------------------------------------------
# Create merged BAM files for replicated samples
#----------------------------------------------

rule merge_replicates:
    output: "02_Bam/{grp}/{sample}-merged.bam"
    input: n_replicates
    shell:
        r"""
        samtools merge - {input} | samtools sort - > {output}
        """

#-----------------------------------------
# Index all BAM files
#-----------------------------------------

rule index_bam:
    input: "02_Bam/{group}/{sample}.bam"
    output: "02_Bam/{group}/{sample}.bam.bai"
    shell: """samtools index {input} {output}"""

#-------------------------------
# Perform Cross-correlation QC
#-------------------------------

# Cross correlation
checkpoint SPP:
    input:
        "02_Bam/{sample}/{sample}{ext}.bam"
    output:
        CC_SCORES_FILE="04_QC/Cross_Correlation/{sample}/{sample}{ext}.cc.qc",
        CC_PLOT_FILE="04_QC/Cross_Correlation/{sample}/{sample}{ext}.cc.plot.pdf"
    threads: 5
    log : "00_Logs/Phantompeakquals/{sample}/{sample}{ext}.log"
    shell:
        r"""
        
        run_spp.R -c={input} -p={threads} -savp={output.CC_PLOT_FILE} -out={output.CC_SCORES_FILE} 2> {log}
        
        sed -r 's/,[^\t]+//g' {output.CC_SCORES_FILE} > temp
        mv temp {output.CC_SCORES_FILE}
        """

# Function to return Fragment length of a SE sample from
# Phantompeakquals output for deepTools

def frag_size(wc):
    # If the sample is PE
    if Files.Strandeness[Files.ID == wc.sample + wc.ext].values[0] == 'PE':
        # determine frag length from mates
        return('')
    # If the sample is SE
    else:
        # Get the output of Phantompeakquals
        checkpoint_output = checkpoints.SPP.get(**wc).output.CC_SCORES_FILE
        df = pd.read_csv(checkpoint_output, sep = '\t', header = None)
        # read in the estimated fragment lenght (3rd collumn)
        est_FragLen = df[2].values[0]
        return(est_FragLen)

#----------------------
# Create genome tracks
#----------------------

rule bam_coverage:
    output: "03_BigWigs/Coverage/{sample}{ext}.bw"
    input: 
        bam = "02_Bam/{sample}/{sample}{ext}.bam",
        bai = "02_Bam/{sample}/{sample}{ext}.bam.bai",
        CC = "04_QC/Cross_Correlation/{sample}/{sample}{ext}.cc.qc"
        
    log: "00_Logs/BamCoverage/{sample}{ext}.log"
    wildcard_constraints: 
        sample = '|'.join(SAMPLES.Name.tolist())
    threads: 10
    params:
        normalize = config['bam_coverage']['normalize'],
        bin = config['bam_coverage']['bin'],
        smooth = config['bam_coverage']['smooth'],
        effectiveGenomeSize = config['effectiveGenomeSize'],
        # Determine fragment length
        extend = frag_size,
        # For pair end reads count only the first mate
        include = lambda wc: '--samFlagInclude 64' if (Files.Strandeness[Files.ID == wc.sample + wc.ext].values[0] == 'PE') else ''
    shell:
        r"""
        bamCoverage -b {input.bam} -o {output} \
        -p {threads} \
        --binSize {params.bin} \
        --smoothLength {params.smooth} \
        --normalizeUsing {params.normalize} \
        --effectiveGenomeSize {params.effectiveGenomeSize} \
        --extendReads {params.extend} {params.include} 2> {log}
        """
# Function to return path to a control of a sample
def macs2_control(wc):
    ID = wc.sample + wc.ext
    
    inputGroup = BAMTARGETS.loc[BAMTARGETS.ID == ID, 'Control_name'].drop_duplicates()
    controlID = BAMTARGETS.loc[BAMTARGETS.ID == ID, 'Control_ID'].drop_duplicates()
        
    return({'control': "02_Bam/" + inputGroup.sum() + '/' + controlID.sum() + ".bam",
            'control_bai': "02_Bam/" + inputGroup.sum() + '/' + controlID.sum() + ".bam.bai"})

#---------------------------------------
# Create Input normalised genome tracks
#---------------------------------------

rule bam_compare:
    output: "03_BigWigs/Input_normalised/{sample}{ext}.InputNormalised.bw"
    input:
        unpack(macs2_control),
        treat = "02_Bam/{sample}/{sample}{ext}.bam",
        treat_bai = "02_Bam/{sample}/{sample}{ext}.bam.bai",
        CC = "04_QC/Cross_Correlation/{sample}/{sample}{ext}.cc.qc"
    log: "00_Logs/BamCompare/{sample}{ext}.log"
    wildcard_constraints: 
        sample = '|'.join(SAMPLES.Name.tolist())
    threads: 10
    params:
        operation = config['bam_compare']['operation'],
        normalize = config['bam_compare']['normalize'],
        bin = config['bam_compare']['bin'],
        smooth = config['bam_compare']['smooth'],
        effectiveGenomeSize = config['effectiveGenomeSize'],
        scaleMethod = config['bam_compare']['scaleFactorsMethod'],
        # Determine fragment length
        extend = frag_size,
        # For pair end reads count only the first mate
        include = lambda wc: '--samFlagInclude 64' if (Files.Strandeness[Files.ID == wc.sample + wc.ext].values[0] == 'PE') else ''
    shell:
        r"""
        bamCompare -b1 {input.treat} -b2 {input.control} -o {output} \
        -p {threads} \
        --operation {params.operation} \
        --binSize {params.bin} \
        --smoothLength {params.smooth} \
        --scaleFactorsMethod {params.scaleMethod} \
        --normalizeUsing {params.normalize} \
        --effectiveGenomeSize {params.effectiveGenomeSize} \
        --extendReads {params.extend} {params.include} 2> {log}
        """

#----------------------------
# Perform QC using deepTools
#----------------------------

# Function to return BigWig files as a input for multiBigwigSummary
def get_group(wc):
    
    # If factor is not all use only samples of one Factor
    if wc.factor != "all":
        grp_indx = Files.ID.str.contains(wc.factor)
        return(list("03_BigWigs/Coverage/" + Files.loc[grp_indx, 'ID'] + ".bw" ))
    else:
        return(list("03_BigWigs/Coverage/" + Files.ID + ".bw" ))

# Create genome-wide matrix per group
rule multiBigWigSummary:
    output: "03_BigWigs/BigWigSummary/{factor}.npz"
    input: get_group
    log: "00_Logs/BigWigSummary/{factor}.log"
    params:
       bins = config['multiBigWigSummary']['bins']
    threads: 10
    shell:
        r"""
        multiBigwigSummary bins -bs {params.bins} -b {input} -o {output}  --smartLabels -p {threads}
        """
# Compute pair-wise correlation per group
rule plot_correlation:
    output:
        mat = "04_QC/Correlation/{factor}.tsv",
        plot = "04_QC/Correlation/{factor}.pdf"
    input: rules.multiBigWigSummary.output
    log: "00_Logs/PlotCorrelation/{factor}.log"
    params:
        method = config['plotCorrelation']['method'],
        colorMap = config['plotCorrelation']['colormap']
    shell:
        r"""
        plotCorrelation -in {input} \
        -c {params.method} \
        -p heatmap \
        -o {output.plot} \
        --skipZeros \
        --removeOutliers \
        --outFileCorMatrix {output.mat} \
        --plotNumbers \
        --colorMap {params.colorMap} 2> {log}
        """

# Function to use different shapes and colors for PCA
# based on factor and treatment

# Compute PCA per group
rule plot_PCA:
    output: 
        plot = "04_QC/PCA/{factor}.pdf",
        data = "04_QC/PCA/{factor}.dat"
    input: rules.multiBigWigSummary.output
    log: "00_Logs/PCA/{factor}.log"
    params:
        top = config['plotPCA']['top']
    shell:
        r"""
        plotPCA -in {input} \
        --plotFile {output.plot} \
        --outFileNameData {output.data} \
        --ntop {params.top} \
        --log2 2> {log}
        """

# Function to determine input bam files for plotFingerprint
def get_bams(wc):
    grp_indx = Files.ID.str.contains(wc.factor)
    bamfiles = list("02_Bam/" + Files.loc[grp_indx, 'Name'] + '/' + Files.loc[grp_indx, 'ID'] + ".bam" )
    baifiles = list("02_Bam/" + Files.loc[grp_indx, 'Name'] + '/' + Files.loc[grp_indx, 'ID'] + ".bam.bai" )
    return({'bam': bamfiles, 'bai': baifiles})

# Plot fingerprint per group
rule plot_fingerprint:
    output: 
        plot = "04_QC/FingerPrint/{factor}.pdf",
        qualmat = "04_QC/FingerPrint/{factor}.qualmat.tab",
        counts = "04_QC/FingerPrint/{factor}.counts.tab"
    input: unpack(get_bams)
    params:
    threads: 10
    run:
        # Determine labels from the paths of bam files
        labels = [i.split('/')[2].replace('.bam', '') for i in input.bam]
        
        # run the shell command
        shell("""
        plotFingerprint -b {input.bam} \
        -plot {output.plot} \
        --labels {labels} \
        --skipZeros \
        -p {threads} \
         --outQualityMetrics {output.qualmat} \
         --outRawCounts {output.counts} \
        """)

#-----------------------------------------------
# Create pseudoreplicates for replicated samples
#-----------------------------------------------

rule create_pseudoreplicates:
    input: "02_Bam/{group}/{sample}.bam"
    output:
        PR1 ="02_Bam/{group}/{sample}_PR1.bam",
        PR2 ="02_Bam/{group}/{sample}_PR2.bam"
    shell:
        r"""
        samtools view -b -s 32.5 {input} > {output.PR1}
        samtools view -b -s 10.5 {input} > {output.PR2}
        """

#------------
# Call peaks
#------------

rule macs2:
    input:
        unpack(macs2_control),
        treatment = "02_Bam/{sample}/{sample}{ext}.bam"
    output: expand("05_Peaks/{{sample}}/{{sample}}{{ext}}{end}", end = ['_peaks.xls', '_peaks.narrowPeak', '_summits.bed' ])
    log : "00_Logs/MACS2/Stringent/{sample}{ext}.log"
    wildcard_constraints: 
        ext = '.+(?<!_PR\d)' # Do not use this method for pseudoreplicates
    params: 
        pair = lambda wc: 'BAMPE' if (design.loc[design.Name == wc.sample, 'Strandeness'].values[0] == 'PE') else 'BAM',
        treshold = config['macs2']['q']
    shell:
        r"""
        macs2 callpeak -t {input.treatment} -c {input.control} \
        --outdir 05_Peaks/{wildcards.sample} \
        -f {params.pair} -g mm -n {wildcards.sample}{wildcards.ext} -q {params.treshold} 2> {log}
        """

#-----------------------------------------------------------
# Calculate Fractions of reads in peaks or blacklist region
#-----------------------------------------------------------

# Calculate Percentage of reads in blacklist regions
rule FRIB:
    output: temp("04_QC/{sample}{ext}_Blacklist.tsv")
    input:
        bam = rules.align.output,
        bai = "02_Bam/{sample}/{sample}{ext}-Raw.bam.bai"
    wildcard_constraints: ext = "-rep\d"
    threads: 3
    shell:
        r"""
        python3 bin/FRIP.py --bam {input.bam} \
        --bed {config[blacklist]} \
        -n {wildcards.sample}{wildcards.ext} \
        -p {threads} > {output}
        """

# Calculate FRIP
rule FRIP:
    output: temp("04_QC/{sample}{ext}_FRIP.tsv")
    input:
        bam = "02_Bam/{sample}/{sample}{ext}.bam",
        bai = "02_Bam/{sample}/{sample}{ext}.bam.bai",
        bed = "05_Peaks/{sample}/{sample}{ext}_peaks.narrowPeak"
    wildcard_constraints: ext = "-rep\d"
    threads: 3
    shell:
        r"""
        python3 bin/FRIP.py --bam {input.bam} \
        --bed {input.bed} \
        -n {wildcards.sample}{wildcards.ext} \
        -p {threads} > {output}
        """

# Concatenate all QC results into single table
rule QC_summary:
    input: expand("04_QC/{sample}_{type}.tsv", sample = design.ID.tolist(), type = ["FRIP", "Blacklist"])
    output: "04_QC/QC_summary.csv"
    run:
        IDs = design.ID.tolist()
        df = pd.DataFrame(columns=("Name", "Blacklist", "FRIP"))
        for ID in IDs:
            FRIP = open("04_QC/" + ID + "_FRIP.tsv", "r").read().strip('\n')
            Blacklist = open("04_QC/" + ID + "_Blacklist.tsv", "r").read().strip('\n')
            
            # Append to df
            df = df.append({'Name': ID, 'Blacklist (%)': Blacklist, 'FRIP (%)': FRIP}, ignore_index=True)
        
        # Save data frame
        df.to_csv(output[0], header = True, index = False)

#-----------------------------------------------
# Call peaks with lower stringency for IDR analysis
#-----------------------------------------------

rule macs2_IDR:
    input:
        unpack(macs2_control),
        treatment = "02_Bam/{sample}/{sample}{ext}.bam"
    output: expand("05_Peaks/{{sample}}/IDR/{{sample}}{{ext}}{end}", end = ['_peaks.xls', '_peaks.narrowPeak', '_summits.bed' ])
    log : "00_Logs/MACS2/IDR/{sample}{ext}.log"
    wildcard_constraints: ext = '-rep\d_PR\d|-merged_PR\d|-rep\d|-merged'
    params: 
        pair = lambda wc: 'BAMPE' if (design.loc[design.Name == wc.sample, 'Strandeness'].values[0] == 'PE') else 'BAM',
        pval = config['macs2_IDR']['p']
    shell:
        r"""
        macs2 callpeak -t {input.treatment} -c {input.control} \
        --outdir 05_Peaks/{wildcards.sample}/IDR \
        -f {params.pair} -g mm -n {wildcards.sample}{wildcards.ext} -p {params.pval} 2> {log}
        """

# Function to return BED files required for IDR analysis
def IDR_input(wc):
    # Pooled pseudoreplicates
    if wc.ext == '_IDR-Np':
        return({'REP1_PEAK_FILE': "05_Peaks/{sample}/IDR/{sample}-merged_PR1_peaks.narrowPeak",
                'REP2_PEAK_FILE': "05_Peaks/{sample}/IDR/{sample}-merged_PR2_peaks.narrowPeak"})
    
    # True replicates
    elif wc.ext == '_IDR-Nt':
        return({'REP1_PEAK_FILE': "05_Peaks/{sample}/IDR/{sample}-rep1_peaks.narrowPeak",
                'REP2_PEAK_FILE': "05_Peaks/{sample}/IDR/{sample}-rep2_peaks.narrowPeak"})
    
    # Pseudoreplicates for rep 1
    elif wc.ext == '_IDR-N1':
        return({'REP1_PEAK_FILE': "05_Peaks/{sample}/IDR/{sample}-rep1_PR1_peaks.narrowPeak",
                'REP2_PEAK_FILE': "05_Peaks/{sample}/IDR/{sample}-rep1_PR2_peaks.narrowPeak"})
    
    # Pseudoreplicates for rep 2
    elif wc.ext == '_IDR-N2':
        return({'REP1_PEAK_FILE': "05_Peaks/{sample}/IDR/{sample}-rep2_PR1_peaks.narrowPeak",
                'REP2_PEAK_FILE': "05_Peaks/{sample}/IDR/{sample}-rep2_PR2_peaks.narrowPeak"})

#-----------------------------------------------
# Irreproducibility Discovery Rate (IDR) analysis
#-----------------------------------------------

# For now the IDR will work only for n=2. 
rule IDR:
    output: 
        out = "05_Peaks/{sample}/IDR/{sample}{ext}_all.narrowPeak",
        REP1_VS_REP2 = "05_Peaks/{sample}/IDR/{sample}{ext}_consistent.narrowPeak"
    input: unpack(IDR_input)
    log: "00_Logs/IDR/{sample}{ext}.log"
    params:
        # Use a tighter threshold for pooled-consistency
        # since pooling and subsampling equalizes the pseudo-replicates
        # in terms of data quality.
        IDR_THRESH = lambda wc: 0.01 if (wc.ext.endswith('Np')) else 0.05
    threads: 5
    run:
        # Check if the overlap between Peaks is > than 20. Otherwise IDR fails.
        import pybedtools
        
        bed1 = pybedtools.BedTool(input.REP1_PEAK_FILE)
        bed2 = pybedtools.BedTool(input.REP2_PEAK_FILE)
        
        if (bed1 + bed2).count() > 20:
            shell(r"""
            idr --samples {input.REP1_PEAK_FILE} {input.REP2_PEAK_FILE} \
            --input-file-type narrowPeak \
            --output-file {output.out} \
            --log-output-file {log} \
            --rank signal.value \
            --soft-idr-threshold {params.IDR_THRESH} --plot --use-best-multisummit-IDR
            
            # Get peaks passing IDR threshold
            
            IDR_THRESH_TRANSFORMED=$(awk -v p={params.IDR_THRESH} 'BEGIN{{print -log(p)/log(10)}}')

            awk 'BEGIN{{OFS="\t"}} $12>='"${{IDR_THRESH_TRANSFORMED}}"' {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' {output.out} \
            | sort | uniq | sort -k7n,7n > {output.REP1_VS_REP2}
            """)
        else:
            shell(r""" 
            touch {output.REP1_VS_REP2}
            touch {output.out}
            echo "Number of overlaping peaks was lower than 20. IDR analysis was not performed" > {log}
            """)

# Calculate Self-Consistency ratio (SCR) and Rescue ratio and return final Peak BED file
rule IDR_QC:
    output: 
        peaks = "05_Peaks/{sample}/IDR/{sample}_finalPeaks.bed",
        QC = temp("04_QC/{sample}_IDR.tsv")
    input:
        Np = "05_Peaks/{sample}/IDR/{sample}_IDR-Np_consistent.narrowPeak",
        Nt = "05_Peaks/{sample}/IDR/{sample}_IDR-Nt_consistent.narrowPeak",
        N1 = "05_Peaks/{sample}/IDR/{sample}_IDR-N1_consistent.narrowPeak",
        N2 = "05_Peaks/{sample}/IDR/{sample}_IDR-N2_consistent.narrowPeak"
    run:
        import pybedtools
        Np = pybedtools.BedTool(input.Np)
        Nt = pybedtools.BedTool(input.Nt)
        N1 = pybedtools.BedTool(input.N1)
        N2 = pybedtools.BedTool(input.N2)
        Np_count = Np.count()
        Nt_count = Nt.count()
        N1_count = N1.count()
        N2_count = N2.count()
        
        # Calculate self-consistency ratio
        if N1_count == 0 or N2_count == 0:
            SCR = 999
        else:
            SCR = max([N1_count, N2_count])/min([N1_count, N2_count])
        
        # Calculate Rescue ratio
        if Np_count == 0 or Nt_count == 0:
            RR = 999
        else:
            RR= max([Np_count, Nt_count])/min([Np_count, Nt_count])
        
        # Determine result
        if SCR < 2 and RR < 2:
            Result = "Ideal"
        elif SCR < 2 or RR < 2:
            Result = "Acceptable"
        elif SCR > 2 and RR > 2:
            Result = "Concerning"
        
        # Export the data to a file
        f = open(output.QC, "w")
        f.write("\t".join([wildcards.sample, str(SCR), str(RR), Result])+"\n")
        f.close()
        
        # Take the longer of Np or Nt as final peak list
        if Np_count > Nt_count:
            Np.saveas(output.peaks)
        else:
            Nt.saveas(output.peaks)

# Condense IDR_QC into single table which will be used for multiQC
rule IDR_QC_summary:
    input: list("04_QC/" + SAMPLES.loc[SAMPLES.N_replicates > 1, 'Name'] + "_IDR.tsv")
    output: "04_QC/IDR_summary.tsv"
    shell:
        r"""
        printf "Sample\tSCR\tRR\tResult\n" | cat - {input} > {output}
        """ 

# Use the narrowPeak file as final peak list for samples with one replicate
rule final_peak_single_rep:
    output: "05_Peaks/Consistent/{sample}_finalPeaks.bed"
    input: lambda wc: "05_Peaks/{sample}/{sample}-rep1_peaks.narrowPeak" if (SAMPLES.loc[SAMPLES.Name == wc.sample,  'N_replicates'].values[0] == 1) else "05_Peaks/{sample}/IDR/{sample}_finalPeaks.bed"
    shell:
        r"""
        ln {input} {output}
        """

#-----------------------------
# Annotate reproducible Peaks 
#----------------------------

# Annotate peaks
# Homer annotation files need to be download first
# perl /home/s1469622/miniconda3/envs/ChIP-seq/share/homer-4.10-0/configureHomer.pl -install mm10

rule annotate_peaks:
    output: "04_QC/Annotation/{sample}_PeakAnnotation.bed"
    input: "05_Peaks/Consistent/{sample}_finalPeaks.bed"
    params:
        genome = config['genome']
    log: "00_Logs/AnnotatePeaks/{sample}.log"
    shell:
        r"""
        annotatePeaks.pl {input} {params.genome} > {output} 2> {log}
        """

#------------------------
# Create annotation plots
#------------------------

# Create plots for annotations
rule plot_anno:
    output:  "04_QC/Annotation/Peak_annotation.summary.txt"
    input: expand("04_QC/Annotation/{sample}_PeakAnnotation.bed", sample = SAMPLES.Name.tolist())
    run:
        Sample_names = ",".join(SAMPLES.Name.tolist())
        inputs = ",".join(input)
        shell(r"""
        ./bin/plot_homer_annotatepeaks.r \
        -i {inputs} -s {Sample_names} \
        -o 04_QC/Annotation -p Peak_annotation
        """)

#--------------------------------------
# Generate summary report using multiQC
#--------------------------------------

rule multiQC:
    input:
            "04_QC/QC_summary.csv",
            "04_QC/IDR_summary.tsv",
            list('04_QC/FingerPrint/' + design.Factor.drop_duplicates() + '.qualmat.tab'),
            list('04_QC/FingerPrint/' + design.Factor.drop_duplicates() + '.counts.tab'),
            list('04_QC/PCA/' + design.Factor.drop_duplicates() + '.dat'),
            list('04_QC/Cross_Correlation/' + Files.Name + '/' + Files.ID + '.cc.qc'),
            list('04_QC/Correlation/' + design.Factor.drop_duplicates() + '.tsv'),
            "04_QC/Correlation/all.tsv",
            "04_QC/PCA/all.dat",
            "04_QC/Annotation/Peak_annotation.summary.txt"
    output: "10_multiQC/multiQC_log.html"
    log: "00_Logs/multiQC.log"
    message: "multiqc for all logs"
    run:
        import yaml
        # Add header to multiqc_config.yaml
        with open('config/multiqc_config.yaml','r') as yamlfile:
            cur_yaml = yaml.safe_load(yamlfile) # Note the safe_load
            cur_yaml['report_header_info'].append(multiqc_header)
        
        if cur_yaml:
            with open('10_multiQC/multiqc_config.yaml','w') as yamlfile:
                yaml.safe_dump(cur_yaml, yamlfile) # Also note the safe_dump
        
        # Run multiQC
        shell("""
        multiqc 00_Logs 01_Fastq/FastQC 04_QC 05_Peaks -o 10_multiQC -c 10_multiQC/multiqc_config.yaml -f -v -n multiQC_log 2> {log}
        """)
