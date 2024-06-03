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


### MicroRNAs
