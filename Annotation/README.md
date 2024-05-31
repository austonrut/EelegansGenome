## This File contains information and step by step instructions for annotating the Edwardsia elegans genome

### Step 1: Repeat Masking and Modeling
Before annotating the genome we need to model the repeats and then mask them. This was done using RepeatModeler and RepeatMasker, respectively. 

For RepeatModler, a new directory was created with only the finished genome fasta file and the bash script to run the following code:

`RepeatModeler  -database edwardsia -recoverDir RM_715597.TueAug11145072023`

RepeatModeler was run once, but I did not give it enough time, to it was run again with the -recoverDir flag to pick up where it left off. This ran for over 100 hours so be sure to give it enought time. 

After RepeatModeler, the output was used to run through RepeatMasker

`RepeatMasker -lib edwardsia-families.fa  Edwardsia.V3P2.final.fasta`

The output from RepeatMasker was then used as input into BRAKER2

### Step 2: BRAKER

`module load braker/2.1.5
module load diamond


braker.pl --cores=29 --species=edwardsia_elegans --genome=Edwardsia.V3P2.final.fasta.masked --prot_seq=prot.metastella.fasta  --bam=../EdV3P2.alignsort.bam \
--prot_seq=/scratch/arutle14/Edwardsia_genome/polishbkr/FinalBraker/RptMdlr/prot.metastella.fasta --etpmode \
--softmasking --AUGUSTUS_CONFIG_PATH=/scratch/arutle14/Edwardsia_genome/polishbkr/config --DIAMOND_PATH=/apps/pkg/diamond/2.0.9/bin \
--ALIGNMENT_TOOL_PATH=${PROTHINT}/bin`
