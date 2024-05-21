# PGGB_starship_pipeline

This pipeline uses currently public and currently in-house assemblies of principally Penicillium and Aspergillus species groups to generate genome graphs  <br />
These genome graphs are then used to extract regions not present in at least 1 other genome <br />
The candidate Regions of Interest (ROIs) are then merged (several steps) and the remaining regions are filtered for a minimum size of 30kb  <br />
All the remaining regions are then annotated  <br />
Following this a set of annotated features found in Starship specific genes are searched for within the proteins in order to extract ROIs with Starship-related genes (SRGs), called here Starship-related regions (SRRs) <br />

The Species groups were selected based on available long-read genomes, genetic distances and interest in relation to domestication (or as a control in the apparent absence of domestication) <br />
Below is an example using a subset of the Fasciculata group in the Penicillium genus containing the species: solitum, crustosum, fuscoglaucum, biforme, caseifulvum and camemberti <br />
This containined a total of 7 genomes (2 for solitum and 1 for each other species) <br />

IMPORTANT <br />
As input the genomes were named in the following format 'species.strain_genomeaccession.fa'. For example 'Psolitum.strain12_ASM1313803v1.fa' or 'Pcaseifulvum.ESE00019.fa' (in the absence of a current public genome accession) <br />
More importantly, as this information is used downstream, the contig names were also modified in a similar format, 'species.strain_contig'. For example '>Psolitum.strain12_JAASRZ010000050' (using the public genome contig/scaffold label) or '>Pcaseifulvum.ESE00019_contig1' <br />


# generate the genome graph and extract the ROIs
```
##to generate the conda/mamba PGGB environment you can run:
#mamba create -n pggb -c bioconda -c conda-forge pggb
##now everything is run in the pggb conda env
conda activate pggb

##First step essentially is to run pggb to generate the graph plus some stats
##here we use a few modifable parameters (some are necessary to modify)
## -n is essentially the number of genomes (will limit the number of paths through a single node to this number therefore attempting to get one path per genome)
## -s is the minimum alignment size so ideally select something over the generic TE sizes, for this example it is ~6-7kb
## -p is the minimum identity...this one is hard to define due to the large distance between genomes...for this analysis is appears safe to go wide of your estimation to be safe
## -m asks for multiqc plots
## -S asks for some wfmash (the aligner) stats
## -t is the number of threads to use during parellelised processes
## -Y is a seperator to generate a sample name per path and this will then be used to avoid self-mappings; so in this case it will be the 'species.strain' name before the underscore '_' that will be given to all contigs from that strain and therefore no within strain mapping will be performed (so we will use the underscore as the seperator)

##just some variables for naming the output files etc
combinedgenome="fasciculata_reduced_less.LR_combined.fa"
combinedgenome2=$( echo $combinedgenome | awk -F ".fa" '{print $1}'  )
##set the number of genomes
genomecount="7"
##set a max divergence beyond what is neccessary
maxdist="75"
threads="30"
soption="30000"
##using default for the kmer size (tried increasing this to 29 as used in yeast examples but for more genome graphs with more diverged species this hindered alignment)
k="19"
##manually selecting and combining the genomes here from a folder called genomes/fasciculata
zcat genomes/fasciculata/Pcaseifulvum.ESE00019.final.renamed.fa.gz genomes/fasciculata/Pfuscoglaucum.ESE00090.final.renamed.fa.gz genomes/fasciculata/Pcamemberti.LCP06093.final.renamed.fa.gz genomes/fasciculata/Pbiforme.LCP05531.final.renamed.fa.gz genomes/fasciculata/Psolitum.strain12_ASM1313803v1.renamed.fa.gz genomes/fasciculata/Psolitum.IBT25940_ASM2882975v1.renamed.fa.gz genomes/fasciculata/Pcrustosum.IBT35664_ASM2882740v1.renamed.fa.gz > ${combinedgenome}
##index the combined genome file
samtools faidx ${combinedgenome}
##nnow build the genome graph
pggb -i ${combinedgenome} -o ${combinedgenome2}.pggb_s${soption}_p${maxdist}_k${k}_G7919_G8069_Y -t ${threads} -p ${maxdist} -s ${soption} -m -S -n ${genomecount} -k ${k} -G 7919,8069 -Y _

##move in the output folder so we can try use this graph to extract Starship regions
cd ${combinedgenome2}.pggb_s${soption}_p${maxdist}_k${k}_G7919_G8069_Y 

##the simpliest path was to use the odgi presence-absence 'PAV' function
##it can produce a matrix with whether regions are present or absent in for each strain
##first we need to generate a file with associates each path (in our cases these are contigs) with a sample (which in our case is always the first string, split by underscores, in the contig name)
odgi paths -i ${combinedgenome}.*.*.*.smooth.final.og -L > ${combinedgenome2}.smooth.final.paths.txt
cut -f 1 -d '_' ${combinedgenome2}.smooth.final.paths.txt > ${combinedgenome2}.smooth.final.samples.txt
paste ${combinedgenome2}.smooth.final.paths.txt ${combinedgenome2}.smooth.final.samples.txt > ${combinedgenome2}.smooth.final.path_and_sample.txt

##Now you can just use bins as done for coverage analysis and then give this to the pav analysis
##here I use a 1kb window
bedtools makewindows -g <(cut -f 1,2 ../${combinedgenome}.fai) -w 1000 > ${combinedgenome2}.w1kb.bed
odgi pav -i ${combinedgenome}.*.*.*.smooth.final.og -b ${combinedgenome2}.w1kb.bed -M -p ${combinedgenome2}.smooth.final.path_and_sample.txt > temp
cat temp | head -n1 > temp2
cat temp | tail -n+2 | sort -k1,1 -k2,2V > temp3
cat temp2 temp3 > ${combinedgenome2}.smooth.final.w1kb.multiple_references.pavs.matrix.tsv
rm temp*
##This system worked much more simply than trying to work with the graph nodes
##then can extract the regions, as done with coverage as well, where we can use the absence of any region in any strain (as this is an all-v-all comp)
##so it melts the matrix
##then WITHIN A STRAIN it combines all 1bp overlapping regions with less than 50% covered, then removes all regions smaller than 10kb (trying to eliminate transposon impact)
##afterwards it does a second merge with 20kb gaps allowed (STILL WITHIN STRAIN); this distance was based on the distribution of the nearest neighbour for regions which indentified a peak around 20kb so is targetting gaps appear more frequent that what we should expect
##these gaps that are merged generally represent transposon impacted regions that break up previously continuous regions
##removing the 10kb regions first though removes the chances of many little fragmnets being joint over large distances during the second merging
##distribution looked at using this : cat ${combinedgenome2}.smooth.final.w1kb.multiple_references.pavs.matrix.tsv | awk '{print $1"\t"$2"\t"$3"\t"$5"\tESE00019\n"$1"\t"$2"\t"$3"\t"$6"\tESE00090\n"$1"\t"$2"\t"$3"\t"$7"\tLCP05531\n"$1"\t"$2"\t"$3"\t"$8"\tLCP06093"}' | awk '{if($4 < 0.5) print}' | bedtools merge -d 1 -c 5 -o distinct -delim ";" | awk '{if($3-$2 > 10000) print}' | awk '{if($1==contig){print $2-pos; pos=$3}else {contig=$1; pos=$3}  }' > ${combinedgenome2}.smooth.final.w1kb.multiple_references.pavs.matrix.dist_neighbour2.txt
##actually perform the extraction and merging to get the Regions Of Interest (ROIs)
for i in $( seq 1 1 ${genomecount} )
do
cat ${combinedgenome2}.smooth.final.w1kb.multiple_references.pavs.matrix.tsv | awk -v row="$i" '{if(NR == 1) {header=$(row+4)} else {print $1"\t"$2"\t"$3"\t"$(row+4)"\t"header}}' | awk '{if($4 < 0.5) print}' | bedtools merge -d 1 -c 5 -o distinct -delim ";" | awk '{if($3-$2 > 10000) print}' | sort -k1,1 -k2,2n  | bedtools merge -d 20000 -c 4 -o distinct -delim ";"
done > ${combinedgenome2}.smooth.final.w1kb.multiple_references.pavs.matrix.ROI.tsv

##Now we can generate a bed file with all these regions and then extract them as done with the coverage analysis
##however we need a cut off point for the size of the regions that will be of real interest
##for this I have now chosen 30kb (in completely contiguous genomes this could probably be increased to ~50kb due to the smallest starships being around this size however 30kb is larger than most regions only filled with transposons and later these regions will be filtered for Starship-related genes anyway)
##after selecting for this cut off we can also merge all the regions into a single call per backbone, keeping track of all the strains within this region and only allowing for a single base gap (previously I allowed for more of a gap but it never made much sense)

##set the minimum ROI size variable
minsizeinsertion="30000"
minsizeinsertion2=$( echo $minsizeinsertion | awk '{print $1/1000}' )
##filter the ROI file and merge accross genomes
cat ${combinedgenome2}.smooth.final.w1kb.multiple_references.pavs.matrix.ROI.tsv | awk -v minsize="$minsizeinsertion" '{if($3-$2 > minsize) print}' | sort -k1,1 -k2,2n | bedtools merge -d 1 -c 4 -o distinct -delim ";"  > ${combinedgenome2}.smooth.final.w1kb.multiple_references.pavs.matrix.ROI.${minsizeinsertion2}kb_min.tsv
##extract the genome regions from these ROIs
cat ${combinedgenome2}.smooth.final.w1kb.multiple_references.pavs.matrix.ROI.${minsizeinsertion2}kb_min.tsv | awk '{print $1":"$2"-"$3}' | while read region
do
samtools faidx ../${combinedgenome} "${region}"
done > ${combinedgenome2}.smooth.final.w1kb.multiple_references.pavs.matrix.ROI.${minsizeinsertion2}kb_min.fa

##in these cases we can certainly be looking at deletions
##however large deletions are rather rare, particularly in coding rich regions
##filtering for some features without knowing the phylogenetic realtionship of the strains could remove important regions. Considering in each graph the phylogenetic relationship was changing, these potential deletions are left in and we can continue to analyse in the absence of potentially biasing the results

##how about some stats about these ROIs, per strain
##here we are interested in a few features
##the number of ROIs (however this can be impacted by genome assembly contiguity)
##the total number of bases in the ROIs (not impacted by contiguity so a better measure for comparing these raw results however again these will be filtered for Starship-related ROIs later)
##we can also look at the average and lax size of the events
echo "strain;count;sum_length;mean_length;max_length" | sed 's/;/\t/g' > ${combinedgenome2}.smooth.final.w1kb.multiple_references.pavs.matrix.ROI.${minsizeinsertion2}kb_min.stats.tsv
cat ${combinedgenome2}.smooth.final.w1kb.multiple_references.pavs.matrix.ROI.${minsizeinsertion2}kb_min.tsv | cut -f1 | awk -F "_" '{print $1}' | sort -u | while read strain
do
grep ^$strain ${combinedgenome2}.smooth.final.w1kb.multiple_references.pavs.matrix.ROI.${minsizeinsertion2}kb_min.tsv | awk -v strain="$strain" '{if($1 ~ strain"_" && max < ($3-$2)) {count++ ; sum=sum+($3-$2); max=($3-$2)} else if($1 ~ strain"_" && max > ($3-$2)) {count++ ; sum=sum+($3-$2)}} END{print strain"\t"count"\t"sum"\t"sum/count"\t"max}'
done >> ${combinedgenome2}.smooth.final.w1kb.multiple_references.pavs.matrix.ROI.${minsizeinsertion2}kb_min.stats.tsv

```
 # Annotate the ROIs with braker2 and funannotate
```

##so now we have our regions of interest across the many strains
##because we do not have annotations for all genomes (depending if public etc) or not uniformly performed annotations at the least; we can now just take our regions of interest and run some annotations on them directly in order to have uniform predictions on all regions (this also saves on annotating entire genomes for each assembly)
##first we need to predict genes
##to do this I used braker installed in a conda env
#mamba create -n braker -c bioconda braker
#conda activate braker
#mamba install -c bioconda -c conda-forge braker2
conda activate braker

##braker requires some other set up steps, here is what I did
###download some augustus config files
#git clone https://github.com/Gaius-Augustus/Augustus
##go to genemark website and download program in order to install it http://exon.gatech.edu/GeneMark/license_download.cgi.
##once downloaded and transfered to current directory
##unpack it
#tar -zxf gmes_linux_64_4.tar.gz
##transer the key to you user home
#cp gmes_linux_64_4/gm_key /home/USER/.gm_key

###we also need some ortholog datasets to feed into braker for training
###this one is the least specific and takes all protein databases for fungi and concatenates them into a large dataset (AFTER SOME BENCHMARKING DONE LATER, FOR THIS DATA, THE EUROTIALES BUSCO DB WOULD PROBABLY HAVE BEEN BETTER, so I recomment doing that instead)
#wget https://v100.orthodb.org/download/odb10_fungi_fasta.tar.gz
#tar -zxf odb10_fungi_fasta.tar.gz
#mv fungi fungi_odb10
#cat fungi_odb10/Rawdata/* > fungi_odb10/refseq_db.fa
#rm -r odb10_fungi_fasta.tar.gz


##now we can run braker annotation
##create some variables to simplify the process of using our ROIs
strain="${combinedgenome2}"
input="${combinedgenome2}.smooth.final.w1kb.multiple_references.pavs.matrix.ROI.${minsizeinsertion2}kb_min.fa"

mkdir ROI_annotation_${minsizeinsertion2}kb_min
cd ROI_annotation_${minsizeinsertion2}kb_min
##reformat the region files, because the names can be too long and some nucleotides (depending on the origin of the genome file) can be softmasked
##make sure all nucleotides are uppercase in the fasta file
##modify the names by just simplifying them to 'strain_regionN'
##also generate a file to keep track of name changes
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}'  ../${input}  | awk -F "." '{if($1 ~ ">") {print ">"$2} else {print}}' | awk -F "_" '{if($1 ~ ">") {print $1}else {print}}' | awk '{if(NR == 1) {strain=$1 ; count=1; print $1":HTR"count} else if($1 ~ ">" && strain==$1){count++ ; print $1":HTR"count} else if($1 ~ ">" && strain!=$1) {$1 ~ ">" && strain=$1 ; count=1; print $1":HTR"count}}' | sed 's/>//g' > temp1
grep '>' ../${input} | sed 's/>//g' > temp2
paste temp2 temp1 > name_mod_association.tsv
rm temp1 temp2

fastafile=$(  echo $input | sed "s/.fa/.UPPERCASE_RENAMED.fa/g" )

awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}'  ../${input}  | awk -F "." '{if($1 ~ ">") {print ">"$2} else {print}}' | awk -F "_" '{if($1 ~ ">") {print $1}else {print}}' | awk '{if(NR == 1) {strain=$1 ; count=1; print $1":HTR"count} else if($1 ~ ">" && strain==$1){count++ ; print $1":HTR"count} else if($1 ~ ">" && strain!=$1) {$1 ~ ">" && strain=$1 ; count=1; print $1":HTR"count} else {print}}' > ${fastafile}
##annotate with braker first to find genes
##need to modify the prot_seq path to the protein database downloaded earlier
##need to modify the AUGUSTUS_CONFIG_PATH and the AUGUSTUS_SCRIPTS_PATH paths to the conda env set up for braker (I was using miniconda3)
##need to modify the GENEMARK_PATH and the PROTHINT_PATH paths to the genemark folder downloaded earlier
##then run
braker.pl --fungus --gff3 --genome=${fastafile} --prot_seq=/PATH/fungi_odb10/refseq_db.fa --epmode --workingdir=braker_output.${strain}.fungi_odb10 --AUGUSTUS_CONFIG_PATH=/PATH/miniconda3/envs/braker/config/ --AUGUSTUS_SCRIPTS_PATH=/PATH/miniconda3/envs/braker/bin/ --GENEMARK_PATH=/PATH/gmes_linux_64_4/ --PROTHINT_PATH=/PATH/gmes_linux_64_4/ProtHint/bin/ --cores=${threads}

##Now we should have predicted genes from braker
##first need to slightly modify the braker output to remove stars placed at the end of the predicted proteins in order to be used by interproscan inside funannotate below
cat braker_output.${strain}.fungi_odb10/augustus.hints.aa | sed 's/\*//g' > braker_output.${strain}.fungi_odb10/augustus.hints.star_mod.aa

##we can actually fully and functionally annotate the genes with protein domains and eggnog annotations using funannotate
##to install funannotate, I did the following
##for this an environement with funannotate was created
#mamba create -n funannotate -c bioconda -c conda-forge funannotate eggnog-mapper
#conda activate funannotate
#run download_eggnog_data.py
##also needed to install a local version of interproscan
#mkdir my_interproscan
#cd my_interproscan
#wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.56-89.0/interproscan-5.56-89.0-64-bit.tar.gz
#wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.56-89.0/interproscan-5.56-89.0-64-bit.tar.gz.md5
#tar -pxvzf interproscan-5.56-89.0-*-bit.tar.gz
#python3 initial_setup.py
#cd ..
##now set up the funnanotate database (this only needs to be done once) using fungi in this case
#funannotate setup -d funannotate_db -b fungi



##activate the env
conda activate funannotate
##set another variable for the threads as the count is set differently in funannotate (it is divisible by eight essentially)
threads3=$( echo $threads | awk '{print $0/8}' | awk -F "." '{print $1}' )
##export the path to the funannotate database as set up before (just using a path to a stored version)
export FUNANNOTATE_DB=/PATH/funannotate_db

##we then run a preliminary step in order to get some interproscan results that aid the annotation step
##for this we use the protein fasta file without the stars at the end 
##need to modify the path to the interproscan downloaded above
funannotate iprscan -i braker_output.${strain}.fungi_odb10/augustus.hints.star_mod.aa -m local -o ${strain}.iprscan -c $threads3 --iprscan_path /PATH/my_interproscan/interproscan-xxx/interproscan.sh

##another preliminary step
##need run egg-nog outside of funannotate due to issue of the path for the eggnogg database (not able to place it in the default position)
##we have already set up the eggnog databases which are downloaded into the non-default position so the path to this db is then given when using emapper.py
##so we just fill in the path to that as the data directory
emapper.py -i braker_output.${strain}.fungi_odb10/augustus.hints.star_mod.aa -o ${strain}.eggnog --data_dir /PATH/eggnog_db --cpu $threads

##now we actually annotate the genes
##so other than using the basic funannotate 'annotate' script, we also feed in the interproscan and eggnog results run previously
funannotate annotate --gff braker_output.${strain}.fungi_odb10/augustus.hints.gff3 --fasta ${fastafile} --species "Penicillium" -o ${strain}.funannotate --strain ${strain} --busco_db fungi --cpus $threads3 --iprscan ${strain}.iprscan --eggnog ${strain}.eggnog.emapper.annotations

##can remove the funannotate and eggnog db as they are large (at the end of everything)
#rm -r funannotate_db
conda deactivate
conda deactivate

```
# Identify the Starship-related genes and extract the ROIs that contain at least one
```
##now identifying the SRGs and their corresponding SRRs

##only take the first transcript annotated
##for this IPR combinations are important and therefore we should looks at combos instead

##a combination of only "IPR009057" && "IPR001005"  or just "IPR001005" and nothing else = the MYB/SANT TF
##previously added the combination of "IPR009057" && "IPR001005" && "IPR017930" (MYB domain) for the MYB/SANT ID but it is too common in other regions of some genomes
##a combination of "IPR013087" (Zinc-finger) && "IPR021842" (DUF3435) or "IPR011010" ("DNA breaking-rejoining enzyme") && "IPR021842" (plus combinations with IPR013762 which is an intergrase-like catalytic domain) and nothing else or "IPR021842" alone = the DUF3435
##any genes with "IPR022198" = the DUF3723
##an interesting case (more common in seems in aspergillus strains) is a combination of the MYB and the IPR013087 (commonly associated with the IPR021842 DUF3435 domain) protein domains = a MYB/SANT-ZnF combo

##the IPRs can be further expanded for other gene identifiers common in starships, starship related genes (SRGs)
##this includes=
## genes containing "ENOG503NY7U" = 'patatin-like phosphatases' as identified in gluck-thaler although this eggnog id is described as 'phosphatidylethanolamine catabolic process'
## genes containing "ENOG503P5Q6" = conidiophore development related genes (CRGs)
## and genes containing "PF20255" = DUF6066 (another DUF with unknown function)

##extract the positions of those genes identified to be SRGs using awk to screen for all the accepted annotations
echo "contig;start;end;gene;Identifying_annotation;SRG_type" | sed 's/;/\t/g' > ${strain}.SRGs_positions.tsv
cat ${strain}.funannotate/annotate_results/Penicillium_${strain}.gff3 | awk -F "\t" '{if($0 ~ "IPR001005" || $0 ~ "IPR021842" || $0 ~ "IPR022198" ) print}' | grep ".t1;" | cut -f9 | awk -F ";" '{print $1}' | sed 's/ID=//g' | sort -u | while read gene
do
IPRs=$( grep $gene ${strain}.funannotate/annotate_results/Penicillium_${strain}.gff3 | cut -f9 | grep InterPro | awk -F "Dbxref" '{print $2}' | awk -F ";" '{print $1}' | tr ',' '\n' | grep InterPro | awk -F ":" '{print $2}' | sort | awk '{sum=sum","$0} END{print sum}' | sed 's/^,//g' )
pos=$( grep $gene ${strain}.funannotate/annotate_results/Penicillium_${strain}.gff3 | grep mRNA | cut  -f1,4-5 | sed 's/_\t/\t/g'  )
echo $pos $gene $IPRs | sed 's/ /\t/g'
done | awk '{if($5 == "IPR001005,IPR009057" || $5 == "IPR001005" || $5 == "IPR021842" || $5 == "IPR013087,IPR021842" || $5 == "IPR011010,IPR021842" || $5 == "IPR011010,IPR013762,IPR021842" || $5 == "IPR011010,IPR013087,IPR013762,IPR021842"   || $5 == "IPR022198" || $5 == "IPR001005,IPR021842" || $5 == "IPR001005,IPR009057,IPR021842" || $5 == "IPR009057,IPR021842" || $5 == "IPR001005,IPR013087" || $5 == "IPR001005,IPR009057,IPR013087" || $5 == "IPR009057,IPR013087" ) {print} else if($5 ~ "IPR021842" && $5 ~ "IPR013087"){print}}' | sort -k1,1 -k2,2n | awk '{if($5 == "IPR001005,IPR009057" || $5 == "IPR001005"  ) {print $0"\tMYB/SANT"} else if($5 == "IPR021842" || $5 == "IPR013087,IPR021842" || $5 == "IPR011010,IPR021842" || $5 == "IPR011010,IPR013762,IPR021842" || $5 == "IPR011010,IPR013087,IPR013762,IPR021842"   ) {print $0"\tDUF3435"} else if($5 ~ "IPR021842" && $5 ~ "IPR013087") {print $0"\tDUF3435"} else if($5 == "IPR022198") {print $0"\tDUF3723"} else if($5 == "IPR001005,IPR021842" || $5 == "IPR001005,IPR009057,IPR021842" || $5 == "IPR009057,IPR021842" ) {print $0"\tDUF3435-MYB/SANT"} else if( $5 == "IPR001005,IPR013087" || $5 == "IPR001005,IPR009057,IPR013087" || $5 == "IPR009057,IPR013087" ) {print $0"\tMYB/SANT"}}' >> ${strain}.SRGs_positions.tsv
cat ${strain}.funannotate/annotate_results/Penicillium_${strain}.gff3 | awk -F "\t" '{if($0 ~ "ENOG503NY7U" ||  $0 ~ "ENOG503P5Q6" ||  $0 ~ "PF20255" ) print}' | grep ".t1;" | cut -f9 | awk -F ";" '{print $1}' | sed 's/ID=//g' | sort -u | while read gene
do
id=$( grep $gene ${strain}.funannotate/annotate_results/Penicillium_${strain}.gff3 | awk '{if($0 ~ "ENOG503NY7U" ){ print "ENOG503NY7U\tPLP"} else if( $0 ~ "ENOG503P5Q6"){ print "ENOG503P5Q6\tCRG"} else if( $0 ~ "PF20255"){ print "PF20255\tDUF6066" } }' )
pos=$( grep $gene ${strain}.funannotate/annotate_results/Penicillium_${strain}.gff3 | grep mRNA | cut -f1,4-5 | sed 's/_\t/\t/g'  )
echo $pos $gene $id | sed 's/ /\t/g'
done >> ${strain}.SRGs_positions.tsv

##we also just need the contig sizes for the genome too 
samtools faidx ${fastafile}
echo "contig;start;end" | sed 's/;/\t/g' > ${strain}.SRGs_positions.contigs_size.tsv
awk '{print $1"\t1\t"$2}' ${fastafile}.fai >> ${strain}.SRGs_positions.contigs_size.tsv

cd ../
##now we can use the annotation to get a better look at the SRRs
##First we can look at gene density as I already observed that fusco regions have very low gene density
echo "assembly;ROI;ROI_gffname;gene_count;size;genes_per_10kb" | tr ';' '\t' > ${combinedgenome2}.smooth.final.w1kb.multiple_references.pavs.matrix.ROI.${minsizeinsertion2}kb_min.gene_density.tsv
cat ROI_annotation_${minsizeinsertion2}kb_min/name_mod_association.tsv | while read regions
do
region1=$( echo $regions | awk '{print $1}' ); region2=$( echo $regions | awk '{print $2}' )
size=$( cat ROI_annotation_${minsizeinsertion2}kb_min/${fastafile}.fai | awk -v region2="$region2" '{if($1 == region2) print $2}'  )
cat ROI_annotation_${minsizeinsertion2}kb_min/${strain}.funannotate/annotate_results/Penicillium_${strain}.gff3 | awk -v region2="$region2" '{if($1 == region2) print}' | grep gene | wc -l | awk -v region1="$region1" -v region2="$region2" -v size="$size" '{print region1"\t"region2"\t"$1"\t"size"\t"($1/(size/10000))}' | awk -F "_" '{print $1"\t"$0}'
done >> ${combinedgenome2}.smooth.final.w1kb.multiple_references.pavs.matrix.ROI.${minsizeinsertion2}kb_min.gene_density.tsv

##now we can extract only the regions that contain SRGs and then recalculate the stats
##extract the list of those with an SRG and then extract those from the generic tsv file with all the ROI info
cat ROI_annotation_${minsizeinsertion2}kb_min/${strain}.SRGs_positions.tsv  | cut -f1 | sort -u | while read region2
do
region1=$( cat ROI_annotation_${minsizeinsertion2}kb_min/name_mod_association.tsv | awk -v region2="$region2" '{if($2 == region2) print $1}' )
cat ${combinedgenome2}.smooth.final.w1kb.multiple_references.pavs.matrix.ROI.${minsizeinsertion2}kb_min.tsv | awk -v region1="$region1" '{if(($1":"$2"-"$3) == region1) print}'
done > ${combinedgenome2}.smooth.final.w1kb.multiple_references.pavs.matrix.ROI.${minsizeinsertion2}kb_min.SRG_related.tsv
##now recalculate the stats and compare it to all the ROIs combined
cat ${combinedgenome2}.smooth.final.w1kb.multiple_references.pavs.matrix.ROI.${minsizeinsertion2}kb_min.stats.tsv | awk '{if(NR==1) {print $0"\tsubset\tcount_proportion\tcount_difference\tsum_length_proportion\tsum_length_difference\tmean_length_difference\tmax_length_difference"}else {print $0"\tall\t1\tNA\t1\tNA\tNA\tNA"}}' > ${combinedgenome2}.smooth.final.w1kb.multiple_references.pavs.matrix.ROI.${minsizeinsertion2}kb_min.SRG_related.stats.tsv
##allow to calculate for strains if they have no SRG-related regions too
cat ${combinedgenome2}.smooth.final.w1kb.multiple_references.pavs.matrix.ROI.${minsizeinsertion2}kb_min.stats.tsv | cut -f1 | tail -n+2 | awk -F "_" '{print $1}' | while read strain
do
grep ^$strain ${combinedgenome2}.smooth.final.w1kb.multiple_references.pavs.matrix.ROI.${minsizeinsertion2}kb_min.SRG_related.tsv | awk -v strain="$strain" '{if($1 ~ strain"_" && max < ($3-$2)) {count++ ; sum=sum+($3-$2); max=($3-$2)} else if($1 ~ strain"_" && max > ($3-$2)) {count++ ; sum=sum+($3-$2)}} END{if(sum > 0){print strain"\t"count"\t"sum"\t"sum/count"\t"max"\tSRG-related"} else {print strain"\t0\t0\t0\t0\tSRG-related"}}'
done | paste - <( grep -v sum_length ${combinedgenome2}.smooth.final.w1kb.multiple_references.pavs.matrix.ROI.${minsizeinsertion2}kb_min.stats.tsv )  | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$2/$8"\t"$2-$8"\t"$3/$9"\t"$3-$9"\t"$4-$10"\t"$5-$11}' >> ${combinedgenome2}.smooth.final.w1kb.multiple_references.pavs.matrix.ROI.${minsizeinsertion2}kb_min.SRG_related.stats.tsv

##extract from the whole fasta file just the regions that are SRG-related
cat ${combinedgenome2}.smooth.final.w1kb.multiple_references.pavs.matrix.ROI.${minsizeinsertion2}kb_min.SRG_related.tsv | awk -v region1="$region1" '{print $1":"$2"-"$3}' | cut -f1 | sort -u | seqkit grep -f - ${combinedgenome2}.smooth.final.w1kb.multiple_references.pavs.matrix.ROI.${minsizeinsertion2}kb_min.fa > ${combinedgenome2}.smooth.final.w1kb.multiple_references.pavs.matrix.ROI.${minsizeinsertion2}kb_min.SRG_related.fa

##can once again calculate gene_density stats but only with the SRG-related tag too
echo "assembly;ROI;ROI_gffname;gene_count;size;genes_per_10kb" | tr ';' '\t' > ${combinedgenome2}.smooth.final.w1kb.multiple_references.pavs.matrix.ROI.${minsizeinsertion2}kb_min.gene_density.SRG_related.tsv
cat ROI_annotation_${minsizeinsertion2}kb_min/${strain}.SRGs_positions.tsv  | cut -f1 | sort -u | while read region2
do region1=$( cat ROI_annotation_${minsizeinsertion2}kb_min/name_mod_association.tsv | awk -v region2="$region2" '{if($2 == region2) print $1}' )
cat ${combinedgenome2}.smooth.final.w1kb.multiple_references.pavs.matrix.ROI.${minsizeinsertion2}kb_min.gene_density.tsv | awk -v region1="$region1" '{if($2 == region1) print}'
done >> ${combinedgenome2}.smooth.final.w1kb.multiple_references.pavs.matrix.ROI.${minsizeinsertion2}kb_min.gene_density.SRG_related.tsv


```

# Combining the results of multiple genome graphs to generalise the impact of domestication


```


```
