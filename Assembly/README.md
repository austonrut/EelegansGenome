## This file contains the step by step and codes used for de-novo assembly of the *Edwardsia elegans* genome from raw nanopore (FAST5) data

### Step 1: Sequenencing
DNA was extracted from a sinlge *E. elegans* individual and sequeneced on a Nanopore minion (LSKxx) and Illumina NextSeq (xxx). 

### Step 2: Basecalling
FAST5 files (both pass & fail) generated from Nanopore sequencing are added to a single directory together (fast5all) and then basecalled using GUPPY.
Basecalled files are put into a new directory (basecallV1.out)

`guppy_basecaller -x "auto" -i ./fast5all/ -s ./basecallV1.out --min_qscore=7 --num_callers 1 -c dna_r10.3_450bps_sup.cfg`

### Step 3: Assembly
Files genereated that pass basecalling are then assembled using wtdbg2

`wtdbg2 -x ont -g 330m -t 16 -i ./basecallV1.out/pass/*fastq -fo EdwRed
wtpoa-cns -t16 -i EdwRed.ctg.lay.gz -fo EdwRed.ctg.fa`

### Step 4: Polish
The assemebled genome from wtdbg2 is then polished with raw illumina reads using Pilon.
To do this, raw Illumina reads are alligned to the assembled genome using bwa, generating a sam file which is then converted to a bam file.
This bam file is used as inputs for Pilon.

`bwa index EdRed3.3polish.fasta`

`bwa mem EdRed3.polish.fasta  E2_S14_L001_R1_001.fastq E2_S14_L001_R2_001.fastq > EdRed3.illalgn.sam`

`samtools view -S -b EdRed3.illalgn.sam  > EdRed3.illalgn.bam`

`samtools sort  EdRed3.illalgn.bam -o EdRed3.illsort.bam`

`samtools index EdRed3.illsort.bam`

`java -Xmx350G -jar ${HOME}/pilon-1.24.jar --genome EdRed3.polish.fasta --bam EdRed3.illsort.bam  --changes  --vcf --diploid --threads 16 --output EdRed3.polish`
