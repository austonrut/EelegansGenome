## This files contains stwp by step how to analysis genomic data

We analysied the genome and protein predictions for orthologous proteins, Macrosynteny, microsynteny, ultra conserved noncoding elements (UCNEs), and microRNAs (miRNAs)

### Orthologous Proteins
Orthofinder was used to identify orthogroups, orthologs, and gene duplications between the newly generated *Edwardsia elegans* genome and the following anthozoans: *Exaiptasia diaphana*, *Actinia tenebrosa*, *Nematostella vectensis*, *Scolanthus callimorphus*, *Acropora millepora*, and *Renilla reniformis*.

`orthofinder -M msa -t 16  -f ./`

#### Single Copy Orthologs
Using Orthofind we can pull of all the shared single copy orthologs between all three Edwardsidae species

`cat Orthogroups.GeneCount.tsv | cut -f1,5,6,8 | awk '$2==$3 && $2==$4 {print $0}'|awk '$2>0' | awk '$2==1 {print $1}'> SingleOgs.EdNvSc.list`

`grep -f 'SingleOGs.EdNvSc.list' Orthogroups.tsv | cut -f1,5,6 > SingleOGs.genenames.EdScUk`

### Macrosynteny


### Microsynteny


### Ultra Conserved Noncoding Elements
Sequences for UCNE were BLAST against the genomes for *E. elegans*, *N. vectensis*, and *S. callimorphus*

`inster blast code here`

Then top hits for each UCNE was used to make files for MacroSyntR

"Bed" files

`cat ucne.Edtest.blast| awk '!a[$1]++' | cut -f1,2,9,10 | awk '{print $2"\t"$3"\t"$4"\t"$1}' > Ed.topucne.bed`
`cat ucne.Nvtest.blast| awk '!a[$1]++' | cut -f1,2,9,10 | awk '{print $2"\t"$3"\t"$4"\t"$1}' > Nv.topucne.bed`
`cat ucne.Sctest.blast| awk '!a[$1]++' | cut -f1,2,9,10 | awk '{print $2"\t"$3"\t"$4"\t"$1}' > Sc.topucne.bed`

"Ortholog" table

`awk '{print $0,$0}' nv.test.list | sed -E 's/(NV.*)( )(NV.*)/\1_ed\t\3_nv/' > EdNV.ucne.tbl`


Only *E. elegans* contigs with more than 1 UCNE were taken for MacroSyntR analysis.



### MicroRNAs
