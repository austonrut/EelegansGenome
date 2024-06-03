## This File contains information and step by step instructions for annotating the Edwardsia elegans genome

To annontate the genome, we first ran BRAKER2 to generate protein predictions, then these predictions were queried using Diamond BLAST to identify the proteins/genes. Protein Family (PFAM) analysis was also done

### Step 1: Repeat Masking and Modeling
Before annotating the genome we need to model the repeats and then mask them. This was done using RepeatModeler and RepeatMasker, respectively. 

For RepeatModler, a new directory was created with only the finished genome fasta file and the bash script to run the following code:

`RepeatModeler  -database edwardsia -recoverDir RM_715597.TueAug11145072023`

RepeatModeler was run once, but I did not give it enough time, to it was run again with the -recoverDir flag to pick up where it left off. This ran for over 100 hours so be sure to give it enought time. 

After RepeatModeler, the output was used to run through RepeatMasker

`RepeatMasker -lib edwardsia-families.fa  Edwardsia.V3P2.final.fasta`

The output from RepeatMasker was then used as input into BRAKER2

### Step 2: BRAKER Protein Predictions
BRAKER2 was run to generate protein predictions. Inputs for predicitions were proteins from *Nematostella vectensis*, the ProtHint metazoan protein database, and raw RNA-seq data from *Edwardsia elegans*

`braker.pl --cores=29 --species=edwardsia_elegans --genome=Edwardsia.V3P2.final.fasta.masked --prot_seq=prot.metastella.fasta  --bam=../EdV3P2.alignsort.bam \
--prot_seq=/scratch/arutle14/Edwardsia_genome/polishbkr/FinalBraker/RptMdlr/prot.metastella.fasta --etpmode \
--softmasking --AUGUSTUS_CONFIG_PATH=/scratch/arutle14/Edwardsia_genome/polishbkr/config --DIAMOND_PATH=/apps/pkg/diamond/2.0.9/bin \
--ALIGNMENT_TOOL_PATH=${PROTHINT}/bin`

### Step 3: BLAST Annotation
First made databases for each Diamond Blast to be run, one for each protein database to be used, NCBI, Uniprot, and *Nematostella vectensis* proteins

`diamond makedb --in allinvert.clean.protein.faa  -d Ncbi
diamond makedb --in uniprot_sprot.fasta -d Unipro
diamond makedb --in NVEC200.20200813.proteins.faa -d JGI`

From there the protein predictions from BRAKER were quired against each data base

`diamond blastp -q rename.augus.hints.aa -d JGI.dmnd -o EdJGI.dblast
diamond blastp -q rename.augus.hints.aa -d Unipro.dmnd -o EdUnipro.dblast
diamond blastp -q rename.augus.hints.aa -d Ncbi.dmnd -o EdNcbi.dblast`

Diamond Blast tables are then sorted by E-values and Percent Similarity, then the top hit for each protein was taken. 

`cat EdUnipro.dblast | sort -k1,1 -k11,11g -k3,3gr | awk '!a[$1]++' > Top.EdUnipro
cat EdJGI.dblast | sort -k1,1 -k11,11g -k3,3gr | awk '!a[$1]++' > Top.EdJGI
cat EdNcbi.dblast | sort -k1,1 -k11,11g -k3,3gr | awk '!a[$1]++' > Top.EdNcbi`

These top hits are then concatenated into a single blast table and the top hit for each protein is taken again, thus giving us a table where each predicted protein is annotated with the best hit based off of three databases

`cat Top.EdUnipro Top.EdJGI Top.EdNcbi | sort -k1,1 -k11,11g -k3,3g | awk '!a[$1]++' > Edele.AnnoBrake.tbl` 

### Step 4: Protein Family Annotation
Protein Families were generated based on the predicted proteins using the program Hmmr. PFams were generated and then filitered to include only hits with an e-value of less than 1. 

`/projects/areitze2_research/PROGRAMS/hmmer-3.1b2/bin/hmmscan \
--cpu 12 --domtblout pfam.domtblout /projects/areitze2_research/PROGRAMS/Pfam-A.hmm rename.augus.hints.aa`

`cat pfam.domtblout | sed -E '1,3d' | awk '$7<1' > pfam.filter.tbl`
