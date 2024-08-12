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
To find macrosytenic regions and generate the oxford dot plot, we first did a recipricol protein blast.

`blastp -query rename.augus.hints.aa  -db uk_jaNemVect1.1_protein.faa -outfmt 6 -out Edp.Nvp.blast`

`blastp -query uk_jaNemVect1.1_protein.faa  -db rename.augus.hints.aa -outfmt 6 -out Nvp.Edp.blast`

From here we take top hits and find common blast hits to get one-to-one protein matches. We do this for each pairwise comparison

`cat Nvp.Edp.blast | awk '$11<1e-05' |awk '!a[$1]++' > Top.NvpEdp.blast`
`cat Top.EdpNvp.blast | cut -f1,2 | sort > Top.EdNv.list`

`comm -1 -2 Top.EdNv.list Top.NvEd.list > EdNv.orthologs.tbl`

Then we use this list to pull out gene locations from a gtf file to generate the input files for MacrosyntR

`cut -f1 EdNv.orthologs.tbl > Ed.ortho2.list`
`cut -f2 EdNv.orthologs.tbl > Nv.ortho2.list`

`grep -f 'Ed.ortho2.list' augustus.hints.gtf | grep -w 'AUGUSTUS.*transcript' | awk '{print $1"\t"$4"\t"$5"\t"$10}' | sed -E 's/(.*)"(.*)";/\1\2/' > Ed2.bed`

`grep -f 'Nv.ortho2.list' uk.prot.chrom | awk '{print $2"\t"$4"\t"$5"\t"$1}' > Nv2.bed`

This leaves us with two files, an orthologs.tbl file which contained one-to-one proteins for two species, and a '.bed' file which isn't really a .bed file but has the name and location of proteins found in the orthologs.tbl file.
We can then use these in the R program MacrosyntR

To compare Macrosynteny between all three Edwardsid species at the same time, we generate a ribbon plot in MacrosyntR using the single copy orthologs generated from Orthofinder
We also `grep` out the locations for these genes the same way as above to generate a .bed file for each species. 

#### MacrosyntR Code
Now that we have generated all the input file we do the rest of the analysis in R, using the following code

`library(dplyr)`
`library(ggplot2)`
`library(egg)`
`library(readr)`
`library(macrosyntR)`

`#For Oxoford Dot Plots`

`my_orth_tbl<- load_orthologs(orthologs_table= "Projects/Data/Edwardsia_Data/GenomeThings/Synteny/macrosyntR/EdSc.orthologs.tbl",
                             sp1_bed="Projects/Data/Edwardsia_Data/GenomeThings/Synteny/macrosyntR/Ed.edsc.mosthits.bed",
                             sp2_bed="Projects/Data/Edwardsia_Data/GenomeThings/Synteny/macrosyntR/Sc.chr.bed")`
                             
`macrosyn_df<- compute_macrosynteny(my_orth_tbl)`

`plot_oxford_grid(my_orth_tbl, sp1_label = "E.elegans", sp2_label = "S.callimorphus", 
                 color_by = "sp1.Chr", dot_size = 1.75, reorder = TRUE, dot_alpha = 0.5)`

`#Generate a list of sytenic regions`
`synt.edsc.list<-get_syntenic_genes(my_orth_tbl)`

`#For Ribbon Plot`

`my_orth_tbl_3<- load_orthologs(orthologs_table= "Projects/Data/Edwardsia_Data/GenomeThings/Synteny/macrosyntR/ScEdNv.3.tbl",
                                bedfiles = c(("Projects/Data/Edwardsia_Data/GenomeThings/Synteny/macrosyntR/Sc3.bed"),
                                              ("Projects/Data/Edwardsia_Data/GenomeThings/Synteny/macrosyntR/Ed3.t25.bed"),
                               ("Projects/Data/Edwardsia_Data/GenomeThings/Synteny/macrosyntR/Nv3.bed")))`

`plot_chord_diagram(my_orth_tbl_3,
                   species_labels = c("S. callimorphus","E. elegans", "N. vectensis"),
                   species_labels_size = 4,
                   color_by = "LGs",
                   ribbons_alpha = 0.5,
                   ribbons_curvature = 0.1,
                   ideogram_height = 5,
                   ideogram_fill = 'white')+
                  theme(legend.position = "none")`

                    
### Microsynteny


### Ultra Conserved Noncoding Elements
Sequences for UCNE (obtained from [https://simrbase.stowers.org/nv_ucne]) were BLAST against the genomes for *E. elegans*, *N. vectensis*, and *S. callimorphus*

`blastn -query nv.ucne.fa  -db ../blast/AnthoBlast/Scal100.genome.fasta  -outfmt 6 -out ucne.Sctest.blast`

Then top hits for each UCNE was used to make files for MacroSyntR

"Bed" files

`cat ucne.Edtest.blast| awk '!a[$1]++' | cut -f1,2,9,10 | awk '{print $2"\t"$3"\t"$4"\t"$1}' > Ed.topucne.bed`
`cat ucne.Nvtest.blast| awk '!a[$1]++' | cut -f1,2,9,10 | awk '{print $2"\t"$3"\t"$4"\t"$1}' > Nv.topucne.bed`
`cat ucne.Sctest.blast| awk '!a[$1]++' | cut -f1,2,9,10 | awk '{print $2"\t"$3"\t"$4"\t"$1}' > Sc.topucne.bed`

"Ortholog" table

`awk '{print $0,$0}' nv.test.list | sed -E 's/(NV.*)( )(NV.*)/\1_ed\t\3_nv/' > EdNV.ucne.tbl`


Only *E. elegans* contigs with more than 1 UCNE were taken for MacroSyntR analysis.
`grep -E 'ctg_0005|ctg_4911|ctg_0813|ctg_0249|ctg_0003|ctg_0099|ctg_0942|ctg_0001|ctg_0101|ctg_0225|ctg_0232|ctg_0388|ctg_0439|ctg_0513|ctg_0520|ctg_0552|ctg_1016|ctg_2070|ctg_3023' Ed.topucne.bed > Ed.top19ucne.bed`

From there we move into R again for MacrosyntR

`my_orth_tbl<- load_orthologs(orthologs_table= "Projects/Data/Edwardsia_Data/GenomeThings/Synteny/macrosyntR/EdNV.ucne.tbl",
                             sp1_bed = "Projects/Data/Edwardsia_Data/GenomeThings/Synteny/macrosyntR/Ed.19ucne.bed",
                             sp2_bed = "Projects/Data/Edwardsia_Data/GenomeThings/Synteny/macrosyntR/Nv.topucne2.bed")`

`macrosyn_df<- compute_macrosynteny(my_orth_tbl)`

`plot_oxford_grid(my_orth_tbl, sp1_label = "E. elegans", sp2_label = "N. vectensis", 
                 color_by = "sp1.Chr", dot_size = 1.75, reorder = TRUE, dot_alpha = 0.5)`


### MicroRNAs
