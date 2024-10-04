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







# Looking for captain-like genes

The captain, defined by the DUF34535/IPR021842 domain has been well described as the conserved element in Starships <br/>
However, what if we wanted to look for other elements, that do not contain a captain, but perhaps a captain-like element, potentially even another Tyrosine recombinase <br/>
In this can we are trying to define not only other captain-like genes, but also other starship-like elements <br/>

First we can evaluate candiates by evaluating all protein domain containing gene in the PAVs <br/>
here there is the assumption that at least one genome in the genome graphs has been fully annotated (XXXXX) in order to evaluate specificity of the domain to the PAV vs the whole genome
```

###a more un-biased means of evalutating the likelihood of certain domains being involved in a HGT like system is to survey all domains
##we can evaluate the domains based on a few stats
##1. their smallest distance to an edge (relative and absolute)
##2. the number of HTR events
##3. the number of genes in total
##4 their frequency in an entire genome (for now we only have one genome annotated that had the most HTRs so we can compare the frequency in this genome in HTRs vs entire)


cd ${dataset}.pggb_s${soption}_p${maxdist}_k${k}_G7919_G8069_Y 
cd ROI_annotation_${minsizeinsertion2}kb_min
##we will first get all interpro domains
##then count the number of HTRs, and genes
##then count how may were present in XXXXX HTRs and how many are present in the whole genome of XXXXX 
##then get their shortest distance to an edge
echo "domain;total_count;shortest_relative_distance_average;shortest_distance_average;HTR_count;gene_count;gene_count_XXXXX;total_XXXXX" | tr ';' '\t' > HTR.protein_domains.stats.tsv
cat ${dataset}.HTRs.gff3 | grep InterPro | awk -F "\t" '{print $9}' | tr ';' '\n' | tr ':' '\n' | tr ',' '\n' | grep ^IPR | sort -u | while read IPR
do
HTRcount=$( grep ${IPR} ${dataset}.HTRs.gff3 | cut -f1 | sort -u | wc -l  )
genecount=$( grep ${IPR} ${dataset}.HTRs.gff3 | awk -F "\t" '{print $9}' | awk -F ";" '{print $1}' | awk -F "=" '{print $2}' | awk -F "." '{print $1}' | sort -u | wc -l )
genecountXX=$( grep ${IPR} ${dataset}.HTRs.gff3 | grep XXXX | awk -F "\t" '{print $9}' | awk -F ";" '{print $1}' | awk -F "=" '{print $2}' | awk -F "." '{print $1}' | sort -u | wc -l  )
totalXX=$( grep ${IPR} ../../../../${dataset}/annotation/3.annotate/3.gff3/XXXX/XXXX.gff3 | awk -F "\t" '{print $9}' | awk -F ";" '{print $1}' | awk -F "=" '{print $2}' | awk -F "." '{print $1}' | sort -u | wc -l )
grep ${IPR} ${dataset}.HTRs.gff3  | while read line
do
HTR=$( echo $line | awk -F " " '{print $1}'  )
pos=$( echo $line | awk -F " " '{print ($4+$5)/2}'  )
cat name_mod_association.tsv | awk -v HTR="$HTR" '{if($2 == HTR) {print $1}}' | tr '-' '\t' | tr ':' '\t' | awk -v pos="$pos" -v line="$line" '{if(($3-$2)-pos >  pos){ print $0"\t"($3-$2)"\t"pos"\t"line} else {print $0"\t"($3-$2)"\t"($3-$2)-pos"\t"line}}'
done | awk -v IPR="$IPR" -v HTRcount="$HTRcount" -v genecount="$genecount" -v genecountXXX="$genecountXXX" -v totalXXX="$totalXXX" '{sum=sum+($5/$4); sum2=sum2+$5; count++} END{print IPR"\t"count"\t"sum/count"\t"sum2/count"\t"HTRcount"\t"genecount"\t"genecount705"\t"total705}'
done >> HTR.protein_domains.stats.tsv

```

Now evaluate these domain containing PAV regions manually using alingments etc with a selected domain (YYYYYY)

```

cd ${dataset}.pggb_s${soption}_p${maxdist}_k${k}_G7919_G8069_Y 

#mamba create -n general -c bioconda -c conda-forge minimap2 samtools bedtools seqkit
conda activate general

domain="YYYYYY"

mkdir HTRs_${domain}
cd HTRs_${domain}

##get list of HTRs with at least one protein with the domain
grep ${domain} ../ROI_annotation_${minsizeinsertion2}kb_min/${dataset}.HTRs.gff3 | cut -f1 | sort -u > HTRs_${domain}.list
##translate the list into the original name
cat HTRs_${domain}.list  | while read HTR; do cat ../ROI_annotation_${minsizeinsertion2}kb_min/name_mod_association.tsv | awk -v HTR="$HTR" '{if($2 ==HTR) print}'; done > HTRs_${domain}.translated.list

##we can extract the proteins and gff3 files for these HTRs
head -n1 ../ROI_annotation_${minsizeinsertion2}kb_min/${dataset}.HTRs.gff3 > HTRs_${domain}.gff3
cat HTRs_${domain}.list | while read HTR
do
cat ../ROI_annotation_${minsizeinsertion2}kb_min/${dataset}.HTRs.gff3 | awk -v HTR="$HTR" '{if($1 == HTR) print}'
done >> HTRs_${domain}.gff3

##now the protein sequences that correspond to these gff3s
cat HTRs_${domain}.gff3 | grep mRNA | awk '{print $9}' | awk -F ";" '{print $1}' | sed 's/ID=//' > HTRs_${domain}.protein.list
seqkit grep -f HTRs_${domain}.protein.list ../ROI_annotation_${minsizeinsertion2}kb_min/${dataset}.HTRs.proteins.fa > HTRs_${domain}.protein.fa

##simplify the gff3 to a bed file with sense, gene name and a label as to whether the gene has any transcript with the domain
echo "contig;start;end;sense;gene;label" | tr ';' '\t' > HTRs_${domain}.genes.bed
cat HTRs_${domain}.protein.list | awk -F "." '{print $1}' | sort -u | while read gene
do
category=$( cat HTRs_${domain}.gff3 | grep "="${gene}"." | awk -v domain="$domain" 'BEGIN{category="NA"} {if($0 ~ domain) {category="found"}} END{if(category=="found") {print "putative-captain"} else {print "NA"}}' )
cat HTRs_${domain}.gff3 | grep "="${gene}"\;" | awk -v gene="$gene" -v category="$category" '{if($3=="gene") {print $1"\t"$4"\t"$5"\t"$7"\t"gene"\t"category}}'
done >> HTRs_${domain}.genes.bed

##now the fasta files
seqkit grep -f HTRs_${domain}.list ../ROI_annotation_${minsizeinsertion2}kb_min/${dataset}.smooth.final.w1kb.multiple_references.pavs.matrix.ROI.20kb_min.UPPERCASE_RENAMED.fa > HTRs_${domain}.fa
##index and create bed file
##index this file
samtools faidx HTRs_${domain}.fa
##then create a bed file
echo "contig;start;end" | tr ';' '\t' > HTRs_${domain}.bed
cut -f1,2 HTRs_${domain}.fa.fai | awk '{print $1"\t1\t"$2}' >> HTRs_${domain}.bed

##align the fasta sequences against itself (but don't allow for self alignment)
minimap2 -cx asm10 --eqx -DP HTRs_${domain}.fa HTRs_${domain}.fa > HTRs_${domain}.minimap_ava.paf



##one way is with the nucleotide sequences, extract the regions of the HTR with a flanking region (in order to show loss of syteny at the edges)
##get the coordinates for each HTR with the protein domain, extend the region by 50kb up and downstream (if possible) and extract from the genome
grep ${domain} ../ROI_annotation_${minsizeinsertion2}kb_min/${dataset}.HTRs.gff3 | cut -f1 | sort -u | while read HTR
do
cat ../ROI_annotation_${minsizeinsertion2}kb_min/name_mod_association.tsv | awk -v HTR="$HTR" '{if($2 == HTR) print}'
done | while read seqs
do
HTR=$( echo $seqs | awk -F " " '{print $2}' )
coord=$( echo $seqs | awk -F " " '{print $1}'  | tr '-' '\t' | tr ':' '\t' | awk '{if(($2-50000) > 0 ) {print $1":"($2-50000)"-"($3+50000)} else {print $1":1-"($3+50000)}}' )
genome=$( echo $seqs | awk -F "_" '{print $1}'  )
samtools faidx ../../genomes/${genome}.fa ${coord} | awk -v HTR="$HTR" '{if($0 ~">") {print ">"HTR} else {print}}'
done > HTRs_${domain}.50kb_plusminus.fa
##index this file
samtools faidx HTRs_${domain}.50kb_plusminus.fa
##then create a bed file
echo "contig;start;end" | tr ';' '\t' > HTRs_${domain}.50kb_plusminus.bed
cut -f1,2 HTRs_${domain}.50kb_plusminus.fa.fai | awk '{print $1"\t1\t"$2}' >> HTRs_${domain}.50kb_plusminus.bed


##we can also look at the genome graphs to get an idea of the actual structure of these HTR regions
##these regions can be complex so have all the alignments visualised at once can help disentangle strata of modification (e.g. clustered insertions)
###considering we are targetting the protein domain ${domain} we can take all region containing at least one gene with the domain
##we can use the file we already created with these locations "HTRs_${domain}.translated.list"
conda activate gg
##used this tutorial https://odgi.readthedocs.io/en/stable/rst/tutorials/extract_selected_loci.html
##need to generate a bed file of the region (+-50kb) of interest wanted to be extracted (alongside all the other regions aligned to it)
##then visualise this region
mkdir HTR_odgi
cat HTRs_${domain}.translated.list | while read line
do
HTR=$( echo $line | awk -F " " '{print $2}' )
coords=$( echo $line | awk -F " " '{print $1}' )
coordsplus=$( echo $coords | tr '-' '\t' | tr ':' '\t' | awk '{if(($2-50000) > 0 ) {print $1":"($2-50000)"-"($3+50000)} else {print $1":1-"($3+50000)}}' )
echo ${coordsplus} | tr ':' '\t' | tr '-' '\t' > temp.bed
##use odgi extract to get the region and all paths associated with the nodes in this region
odgi extract -i ../${dataset}.combined.fa.*.*.*.smooth.final.og -o HTR_odgi/${HTR}.og -b temp.bed -c 0 --threads ${threads} -P
##now visualise after re-sorting
## -s is the character used, to group all contig names before that character as the same genome and colour them differently (will highlight if any genomes have multiple alignments here)
odgi sort -i HTR_odgi/${HTR}.og -o - -O | odgi viz -i - -o HTR_odgi/${HTR}.png -s '_' -P
done
rm temp.bed


##in some cases there are clean PAVs and some other cases are more complicated
##for exmaple XXXXX:53
##it looks like there is two close by insertions (split by a small common region) that were merged into one call. This mimics what is seen in plots of the region and it's ${domain} containing genes
##we can try plot this now using the pggb alignment paf file and isolating these HTRs, basically to produce the odgi image but with the ability to add annotations too
##first generate a bed file of all the genomes and their contigs
echo "contig;start;end" | tr ';' '\t' > ../../${dataset}.combined.bed
cut -f1,2 ../../${dataset}.combined.fa.fai | awk '{print $1"\t1\t"$2}' >> ../../${dataset}.combined.bed
##adjust all gene positions to the whole genome
echo "contig;start;end;sense;gene;label" | tr ';' '\t' > HTRs_${domain}.genes.genome_positions.bed
tail -n+2 HTRs_${domain}.genes.bed | while read line
do
HTR=$( echo $line | awk  -F " " '{print $1}' )
HTR2=$( cat HTRs_${domain}.translated.list | awk -v HTR="$HTR" '{if($2 == HTR) print $1}' )
HTR2contig=$( echo $HTR2 | awk -F ":" '{print $1}' )
HTR2start=$( echo $HTR2 | awk -F ":" '{print $2}' | awk -F "-" '{print $1}' )
echo "" | awk -v line="$line" '{print line}' | awk -v HTR2contig="$HTR2contig" -v HTR2start="$HTR2start" '{print HTR2contig"\t"($2+HTR2start)"\t"($3+HTR2start)"\t"$4"\t"$5"\t"$6}'
done >> HTRs_${domain}.genes.genome_positions.bed


##the pggb alignment is not very clear based on initial observation (they appear too filtered and aggregated to see refined differences (probably a fault of trying to simplify the graph)) so we can try an alternative route
##extract the regions associated with each HTR, as found by odgi, and realign these regions trying to be more refined
cat HTRs_${domain}.translated.list | while read line
do
HTR=$( echo $line | awk -F " " '{print $2}' )
odgi paths -L -i HTR_odgi/${HTR}.og > HTR_odgi/${HTR}.odgi_paths.txt
cat HTR_odgi/${HTR}.odgi_paths.txt | while read region
do
samtools faidx ../../${dataset}.combined.fa "${region}" 
done > HTR_odgi/${HTR}.odgi_paths.fa
samtools faidx HTR_odgi/${HTR}.odgi_paths.fa
echo "contig;start;end" | tr ';' '\t' > HTR_odgi/${HTR}.odgi_paths.bed
cut -f1,2 HTR_odgi/${HTR}.odgi_paths.fa.fai | awk '{print $1"\t1\t"$2}' >> HTR_odgi/${HTR}.odgi_paths.bed
minimap2 -cx asm5 --eqx -DP HTR_odgi/${HTR}.odgi_paths.fa HTR_odgi/${HTR}.odgi_paths.fa > HTR_odgi/${HTR}.odgi_paths.minimap_ava.paf
done

##as above but using nucmer instead, then converting the delta file to a paf file in order to be visualised with gggenomes
##nucmer appears to work better when trying to define these PAVs, otherwise the minimap chaining aligning algorithm appears to force poor long alignemnts through the gap
#mamba create -n mummer4 bioconda::mummer4
conda activate mummer4

cat HTRs_${domain}.translated.list | while read line
do
HTR=$( echo $line | awk -F " " '{print $2}' )
nucmer --maxmatch --delta HTR_odgi/${HTR}.odgi_paths.nucmer.delta --nosimplify HTR_odgi/${HTR}.odgi_paths.fa HTR_odgi/${HTR}.odgi_paths.fa
conda activate general
paftools.js delta2paf HTR_odgi/${HTR}.odgi_paths.nucmer.delta > HTR_odgi/${HTR}.odgi_paths.nucmer.paf
conda deactivate
done

conda deactivate


##adjust coords of genes to these regions (+-50kb edges)
cat HTRs_${domain}.translated.list | while read line
do
HTR=$( echo $line | awk -F " " '{print $2}' )
echo "contig;start;end;sense;gene;label" | tr ';' '\t' > HTR_odgi/${HTR}.odgi_paths.genes_position.tsv
contig=$( echo $line | awk -F " " '{print $1}' | awk -F ":" '{print $1}' )
HTR2mod=$( cat HTR_odgi/${HTR}.odgi_paths.txt | awk -v contig="$contig" '{if($0 ~ contig) print}' )
tail -n+2 HTRs_${domain}.genes.bed | awk -v HTR="$HTR" -v HTR2mod="$HTR2mod" '{if($1 == HTR) print HTR2mod"\t"($2+50000)"\t"($3+50000)"\t"$4"\t"$5"\t"$6}' >> HTR_odgi/${HTR}.odgi_paths.genes_position.tsv
done 

```


# Now use R to plot these data

```{r setup, include=FALSE}
library(ggExtra)
library(ggpubr)
library(gggenomes)
```

## first getting the candiate protein domain

```{r}

##trying to define whether a gene is specific to a HTR region is difficult by eye
##can look at all domains and whether they have features indicative of being at a starship edge and being specific to HTR regions

##read in file containing some stats
IPRtest=read.csv(file="pathto/${dataset}.pggb_s${soption}_p${maxdist}_k${k}_G7919_G8069_Y/ROI_annotation_${minsizeinsertion2}kb_min/HTR.protein_domains.stats.tsv", sep='\t', header=T)

##here we plot a few important stats and highight those that may be of importance
##this includes:
##xaxis = gene_count_XXXXX/total_XXXXX
##          This measure tells us the number of genes with the protein domain of interest were found in the HTRs versus the entire genome for XXXXX. We expect that the majority of the genes whould be HTR specific therefore a ratio of 1 is ideal. However we can assume that not all HTRs are detected and that for now proteins with the same domain but very different may be present.
##yaxis = shortest_relative_distance_average
##          This tells us how close on average is the gene to an edge, relative to the entire size of the HTR. We expect a captain like gene to be closer to an edge compared to other genes.
##size = gene_count/HTR_count
##          This is a measure of the number of genes with the protein domain versus the number of HTRs that at least one gene is found it. It tells us about frequency of the gene in each HTR, ideally it would be a clean 1:1 but this is not always the case even in clean Starships.

##we also filter out domains that:
##    appear in <=2 genes
##    appear in <=2 HTRs
##    appear <=2 times in the whole genome
##    have a higher number of genes detected in the HTR annotation vs whole genome (repeat elements are not annotated in the whole genome dataset)
##    have a ratio of gene_count:HTR_count >=3 (i.e. many genes in few HTRs)
##this cleans up a lot of spurious, uninformative domains

IPRtest2=subset(IPRtest, total_XXXXX > 2 & gene_count_XXXXX > 2 & gene_count > 2 & HTR_count > 2 & gene_count/HTR_count < 3 & gene_count_XXXXX/total_XXXXX <=1)

plot1=ggplot()+
  geom_point(data=IPRtest2, aes(x=(gene_count_XXXXX/total_XXXXX), y=shortest_relative_distance_average, size=gene_count/HTR_count), alpha=0.5)+
  xlim(0,1)+
  geom_point(data=subset(IPRtest2 , gene_count_XXXXX/total_XXXXX >0.5 & shortest_relative_distance_average < 0.25), aes(x=(gene_count_XXXXX/total_XXXXX), y=shortest_relative_distance_average, size=gene_count/HTR_count), colour="red", show.legend = F)+
  geom_density_2d(data=IPRtest2, aes(x=(gene_count_XXXXX/total_XXXXX), y=shortest_relative_distance_average), show.legend = F)+theme_pubr()+geom_hline(linetype="dotted", colour="red", yintercept = 0.25)+geom_vline(linetype="dotted", colour="red", xintercept = 0.5)+xlab("Total proportion in PAVs (XXXXX)")+ylab("Average shortest relative distance to a PAV edge")+ylim(0,0.5)
ggMarginal(plot1, type = "histogram", bins = 100, fill = "steelblue")

```

Now plot PAVs containing candidate domains using gggenomes

```{r, out.width='100%'}
domain="YYYYYY"

##read in bed file of HTRs
HTRs=read.csv("pathto/${dataset}.pggb_s${soption}_p${maxdist}_k${k}_G7919_G8069_Y//HTRs_${domain}/HTRs_${domain}.bed", sep='\t', header=T)
#HTRs=read_seqs("pathto/${dataset}.pggb_s${soption}_p${maxdist}_k${k}_G7919_G8069_Y/HTRs_${domain}/HTRs_${domain}.fa")
##make another header, copy of contig, but called seq_id
HTRs$seq_id = HTRs$contig
##also length using the end position (all HTRs are relative to 0)
HTRs$length = HTRs$end


##plot just the contigs without any features
gggenomes(seqs=HTRs) +geom_seq()+geom_seq_label()

genes=read.csv("pathto/${dataset}.pggb_s${soption}_p${maxdist}_k${k}_G7919_G8069_Y//HTRs_${domain}/HTRs_${domain}.genes.bed", header=T, sep='\t')
genes$seq_id = genes$contig
genes$length= genes$end-genes$start
genes$strand=genes$sense

##plot sequences with genes on top
##wrap the facet so each has its own x-axis scale (only work around with gggenomes I could find), making it easier to see the smaller elements and the putative captain
gggenomes(genes=genes, seqs=HTRs)+
  geom_seq()+geom_gene(aes(fill=label), stroke=0)+
  geom_seq_label()+
  scale_fill_manual(values = c("red","grey"))+
  theme(legend.position="top", legend.box = "horizontal")+
  facet_wrap(seq_id ~ ., scales="free")


##read in the paf all vs all alignment
links=read_links("pathto/${dataset}.pggb_s${soption}_p${maxdist}_k${k}_G7919_G8069_Y/HTRs_${domain}/HTRs_${domain}.nucmer.paf")

##plot the HTRs, the genes and the nucleotides alignments from nucmer
##HTRs are manually rearranged here so similar HTRs (or identical) are side by side
gggenomes(genes=genes, seqs=HTRs, links=subset(links, map_length > 10000 & seq_id != seq_id2), adjacent_only = F) %>%
  pick() %>%
  gggenomes::sync() %>%
  gggenomes::flip(3)+
  geom_link(size=0.1, colour="black", alpha=0.25)+
  geom_seq()+geom_gene(aes(fill=label), stroke=0)+
  geom_seq_label()+
  scale_fill_manual(values = c("red","grey"))+
  theme(legend.position="bottom", legend.box = "horizontal")


```

## Refined alignments of the regions <br/>
As the genome graph alignments were designed to give a good global view of synteny, the refined edges and smaller alignments were lost <br/>
therefore we can then take the regions found as similar in the genome graph and realign with nucmer <<br/>
we can then plost this instead

##example with HTR:ZZ

```{r, out.width='100%'}

###a better way is using the regions to redo the alignment in order to be more accurate with these local regions (can be less stringent than in the genome graph due to the amount of similarity of disparate regions)
##below is an example for HTR:ZZ 

##get the regions of the HTR (+-50kb) bed file
contigs=read.csv("pathto/${dataset}.pggb_s${soption}_p${maxdist}_k${k}_G7919_G8069_Y/HTRs_${domain}/HTR_odgi/HTR:ZZ.odgi_paths.bed", sep='\t', header=T)
contigs$seq_id = contigs$contig
contigs$length = contigs$end


##get the minimap2 paf file for the links
links=read_links("pathto/${dataset}.pggb_s${soption}_p${maxdist}_k${k}_G7919_G8069_Y/HTRs_${domain}/HTR_odgi/HTR:ZZ.odgi_paths.nucmer.paf")

##get the gene positions (with the IPR052925 labelled as putative-captain)
genes=read.csv("pathto/${dataset}.pggb_s${soption}_p${maxdist}_k${k}_G7919_G8069_Y/HTRs_${domain}/HTR_odgi/HTR:ZZ.odgi_paths.genes_position.tsv", sep='\t' , header=T)
genes$seq_id = genes$contig
genes$length= genes$end-genes$start
genes$strand=genes$sense


##plot it...
##may need to rearrange the order of the regions manually
gggenomes(seqs=contigs, links=subset(links), genes=genes)%>%
    pick() %>%
    sync() %>%
    flip()+
    geom_link(aes(fill=(map_match/map_length)), offset = c(0,0), size=0.1)+
    geom_seq()+
    geom_seq_label()+
    geom_gene(aes(colour=label), fill="white")+
    scale_color_manual(values="red","grey")+
    theme(legend.position = "top")


```



