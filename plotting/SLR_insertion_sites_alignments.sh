conda activate general
cd SLR_plots
mkdir between_SLRs_plus

##set the flank for alignment
flank=50000
flank2=$( echo $flank | awk '{print $0/1000}' )
##get all HTRs with an exact match in another HTR and create a txt file of this cluster
cat ../all_SRG_related_ROIs.renamed.NR/EXACTMATCHES.${identity}_${overlap}.tsv ../all_SRG_related_ROIs.renamed.NR/OFINTEREST.${identity}_${overlap}.tsv | cut -f4 | sort -u | while read HTR3
do
HTR=$( cat ../all_SRG_related_ROIs.name_mod_association.tsv | awk -v HTR3="$HTR3" '{if($2 == HTR3) print $1}' )
contig=$( echo $HTR | awk -F ":" '{print $1}' )
start=$( echo $HTR | awk -F ":" '{print $2}' | awk -F "-" '{print $1}' )
end=$( echo $HTR | awk -F ":" '{print $2}' | awk -F "-" '{print $2}' )
HTR2=$( cat  ../*.LR_combined.pggb_s*_p*_k19_G7919_G8069_Y/ROI_annotation_30kb_min/name_mod_association.tsv | agrep -2 "${HTR}" | awk -v contig="$contig" '{if($1 ~ contig":") print $2}' )
echo "${HTR2}" > between_SLRs_plus/cluster_${HTR2}.txt
##get all HTRs that were exactly the same
cat ../all_SRG_related_ROIs.renamed.NR/EXACTMATCHES.${identity}_${overlap}.tsv ../all_SRG_related_ROIs.renamed.NR/OFINTEREST.${identity}_${overlap}.tsv  | awk -v HTR3="$HTR3" '{if($4 == HTR3 && $7 != "") print $7}' | sort -u | while read HTR4
do
HTR5=$( cat ../all_SRG_related_ROIs.name_mod_association.tsv | awk -v HTR4="$HTR4" '{if($2 == HTR4) print $1}' )
contig2=$( echo $HTR5 | awk -F ":" '{print $1}' )
##search for the corresponding HTR (allowing for a single missmatch using agrep, this is becuase at one point I reran some genomegraphs and didn't rerun this part and sometimes the differences is 1kb and would rather more matches than less...can be ignored )
HTR6=$( cat  ../*.LR_combined.pggb_s*_p*_k19_G7919_G8069_Y/ROI_annotation_30kb_min/name_mod_association.tsv | agrep -2 "${HTR5}"  | awk -v contig2="$contig2" '{if($1 ~ contig2":") print $2}' )
echo "${HTR6}" >> between_SLRs_plus/cluster_${HTR2}.txt
done
##get other genomes that apparently don't have the region
others=$( cat  ../*.LR_combined.pggb_s*_p*_k19_G7919_G8069_Y/*.pavs.matrix.ROI.*kb_min.tsv | awk -v contig="$contig" -v start="$start" -v end="$end" '{if($1 == contig && $2 == start || $1 == contig && $3 == end) print $4}' | tr ';' '\n'  )
echo "${others}" > between_SLRs_plus/cluster_${HTR2}.absent.txt
done
rm between_SLRs_plus/cluster_.txt

##loop through the HTRs in each cluster text file and extract the information for the alignment etc
ls between_SLRs_plus/cluster_*.txt | grep -v absent | while read clusterfile
do
cluster=$(  echo "${clusterfile}" | awk -F "/" '{print $NF}' | sed 's/.txt//' )
rm between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.fa
echo "contig;start;end;sense;gene;label" | tr ';' '\t' > between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.genes.bed
cat "${clusterfile}" | while read HTR2
do
##get the alternate HTR name
HTR=$( cat  ../*.LR_combined.pggb_s*_p*_k19_G7919_G8069_Y/ROI_annotation_30kb_min/name_mod_association.tsv | awk -v HTR2="$HTR2" '{if($2 == HTR2) print $1}' )
genome=$( echo "${HTR}" | awk -F "_" '{print $1}' )
contig=$( echo "${HTR}" | awk -F ":" '{print $1}' )
start=$( echo "${HTR}" | awk -F ":" '{print $2}' | awk -F "-" '{print $1}' )
end=$( echo "${HTR}" | awk -F ":" '{print $2}' | awk -F "-" '{print $2}' )
startmod=$( echo "${start}" | awk -v flank="$flank" '{print $0-flank}' )
endmod=$( echo "${end}" | awk -v flank="$flank" '{print $0+flank}' )
##need the get the SLR annotations in order to add them to the plot
##need to modify the positions by adding the flank to each start and end
##need to also highlight the SRGs
cat ../*.LR_combined.pggb_s*_p*_k19_G7919_G8069_Y/ROI_annotation_30kb_min/*/annotate_results/*.gff3 | grep gene | awk -v HTR2="$HTR2" -v flank="$flank" -v start="$start" '{if($1 == HTR2 && start < flank) {print $1"\t"$4+start"\t"$5+start"\t"$7"\t"$9} else if ($1 == HTR2 && start > flank) {print $1"\t"$4+flank"\t"$5+flank"\t"$7"\t"$9}}' | while read line2
do
contig2=$( echo "${line2}" | awk '{print $1}' )
start2=$( echo "${line2}" | awk '{print $2}' )
end2=$( echo "${line2}" | awk '{print $3}' )
sense=$( echo "${line2}" | awk '{print $4}' )
geneid=$( echo "${line2}" | awk -F "\t" '{print $5}' | awk -F ";" '{print $1}' | awk -F "=" '{print $2}' )
label=$( cat ../*.LR_combined.pggb_s*_p*_k19_G7919_G8069_Y/ROI_annotation_30kb_min/*.LR_combined.SRGs_positions.tsv | awk -v HTR2="$HTR2" -v geneid="$geneid" 'BEGIN{label = "NA"} {if($1 == HTR2 && $4 ~ geneid".t" ) {label = $6} } END{print label}' )
##print the output
echo "${contig2};${start2};${end2};${sense};${geneid};${label}" | tr ';' '\t' >> between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.genes.bed
done
##add the captain genes specifically added by starfish (unique to it)
cat ../*.LR_combined.pggb_s*_p*_k19_G7919_G8069_Y/starfish.captains.starfish_unique.bed | grep ^"${contig}" | awk -v contig="$contig" -v start="$start" -v end="$end" '{ if($1 == contig && $2 == start && $3 == end) print}' | while read line
do
cat ../*.LR_combined.pggb_s*_p*_k19_G7919_G8069_Y/starfish_results/geneFinder/*.bed | grep "${line}" | awk -v contig="$contig" -v start="$start" -v flank="$flank" -v HTR2="$HTR2" '{if( $1 == contig && $2 > start)print start"\t"HTR2"\t"$2-start+flank"\t"$3-start+flank"\t"$6"\t"$4"\tDUF3435"}' 
done  >> between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.genes.bed


##get the genome and extract the contig
zcat ../genomes/*/${genome}_*.fa.gz ../genomes/*/${genome}.final.renamed.fa.gz ../genomes/*/${genome}.renamed.fa.gz > ${genome}.fa
samtools faidx ./${genome}.fa ${contig}:${startmod}-${endmod} | awk -v HTR2="$HTR2" '{if($0 ~ ">" ) {print ">"HTR2} else {print}}' >> between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.fa
rm ./${genome}.*

done


##create header for bed file of regions extracted by starship etc, i.e. the aligned regions
echo "contig;start;end;fullname" | tr ';' '\t' > between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.bed



## align genomes without the insertion to the HGT then extract all that align more than 30kb (one full flank) hoping to find the added flanks of the starship
## then take that contig, take all alignments greater than 20kb and extrac the max and min of those alignments, hopefully surrounding the insertion
##then extract just this part of the contig and add to the fasta file for all_v_all alignment
HTR2=$( head -n1 between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.fa | sed 's/>//g' )
samtools faidx between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.fa "${HTR2}" > between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.temp.fa
rm between_SLRs_plus/${cluster}.othergenomes.fa
rm between_SLRs_plus/${cluster}.othergenomes.contigs.fa
cat between_SLRs_plus/${cluster}.absent.txt | while read othergenome
do
zcat ../genomes/*/${othergenome}_*.fa.gz ../genomes/*/${othergenome}.final.renamed.fa.gz ../genomes/*/${othergenome}.renamed.fa.gz >> between_SLRs_plus/${cluster}.othergenomes.fa
done
nucmer --maxmatch --minmatch 100 --delta  between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.othergenomes.nucmer.delta between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.temp.fa between_SLRs_plus/${cluster}.othergenomes.fa
#nucmer --maxmatch --minmatch 100 --delta  between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.othergenomes.nucmer.delta between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.fa between_SLRs_plus/${cluster}.othergenomes.fa
conda activate general
paftools.js delta2paf between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.othergenomes.nucmer.delta > between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.othergenomes.nucmer.paf
cat between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.othergenomes.nucmer.paf | awk '{if($11 > 30000 ) print}' | cut -f1 | sort -u | while read tempcontig
do
cat between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.othergenomes.nucmer.paf | awk -v tempcontig="$tempcontig" '{if($1 == tempcontig) sum=sum+$11} END{print tempcontig"\t"sum}'
done | sort -k2nr | head -n2 | awk '{print $1}' | while read insertioncontig
do
strain=$( echo "${insertioncontig}" | awk -F "_" '{print $1}' )
edges=$( cat between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.othergenomes.nucmer.paf | awk '{if($11 > 20000) print}' | awk -v insertioncontig="$insertioncontig" '{if($1==insertioncontig) print}' | awk -v flank="$flank" 'BEGIN{max=0; min=99999999999999} {if($4 > max) {max=$4}; if($3 < min) {min=$3}} END{print min-flank"\t"max+flank}' | awk '{if($1 < 0) {print 0"\t"$2} else {print}}' )
edges2=$( echo "${edges}" | awk '{print $1"-"$2}' )
size=$( echo "${edges}" | awk '{print $2-$1}' )
samtools faidx between_SLRs_plus/${cluster}.othergenomes.fa "${insertioncontig}:${edges2}" >> between_SLRs_plus/${cluster}.othergenomes.contigs.fa
echo "${insertioncontig}:${edges2};1;${size};NA" | tr ';' '\t' >> between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.bed
done 
rm between_SLRs_plus/${cluster}.othergenomes.fa 
rm between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.temp.fa



##add to the bed file of the SLRs that will be aligned
samtools faidx between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.fa
cut -f1,2 between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.fa.fai | awk '{print $1"\t1\t"$2}' | while read something
do
name=$( echo "${something}" | awk '{print $1}' )
name2=$( cat  ../*.LR_combined.pggb_s*_p*_k19_G7919_G8069_Y/ROI_annotation_30kb_min/name_mod_association.tsv | awk -v name="$name" '{if($2 == name) print $1}'  )
echo "${something};${name2}" | tr ';' '\t'
done >> between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.bed

##create bed file of the region of actual SLR (little more complicated incase there is not enough space for the flanking regions, therefore need to figure that out and adjust for it)
echo "contig;start;end" | tr ';' '\t' > between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.SRG.bed
cut -f1 between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.fa.fai | while read HTRagain
do
start=$( cat  ../*.LR_combined.pggb_s*_p*_k19_G7919_G8069_Y/ROI_annotation_30kb_min/name_mod_association.tsv  | grep $HTRagain | awk -v HTRagain="$HTRagain" '{if($2 == HTRagain)print $1}' | awk -F ":" '{print $2}' | awk -F "-" '{print $1}' ) 
size=$( cat  ../*.LR_combined.pggb_s*_p*_k19_G7919_G8069_Y/ROI_annotation_30kb_min/name_mod_association.tsv  | grep $HTRagain | awk -v HTRagain="$HTRagain" '{if($2 == HTRagain)print $1}' | awk -F ":" '{print $2}' | awk -F "-" '{print $2-$1}' ) 
grep ^${HTRagain} between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.fa.fai | awk -v HTRagain="$HTRagain" '{if($1 == HTRagain) print}' | awk -v start="$start" -v flank="$flank" -v size="$size" '{if(start < flank) {print $1"\t"start"\t"start+size} else {print $1"\t"flank"\t"flank+size}}' 
done  >> between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.SRG.bed


##add the contigs missing the elements to the file before alignment
cat between_SLRs_plus/${cluster}.othergenomes.contigs.fa >> between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.fa

##generate all vs all alignments for each set
nucmer --maxmatch --minmatch 100 --delta  between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.nucmer.delta between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.fa between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.fa
conda activate general
paftools.js delta2paf between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.nucmer.delta > between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.nucmer.paf

## change name from  horizontally transferred regions (HTR) to Starship-like regions (SLR)
sed -i 's/:HTR/:SLR/g' between_SLRs_plus/${cluster}*

##automate the production of an R script using gggenomes to plot the alignment
##then use gggenome with R script to create the plots
cat ../gggenomes_skeleton.R | sed "s/GENOME/${cluster}/g" | sed "s/FLANK2/${flank2}/g" | sed "s/within_genomes/between_SLRs_plus/g" > between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.R
Rscript between_SLRs_plus/${cluster}.SRGs.${flank2}kb_flank.R

done 
