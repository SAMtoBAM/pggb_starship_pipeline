overlap="85"
identity="85"
gap="15000"

#generate a file output region which be be filled with intermediate files such as the BLAST alignments etc
mkdir all_SRG_related_ROIs.renamed.NR
#remove the summary output file incase the same output directory is being used again
rm all_SRG_related_ROIs.renamed.NR/OFINTEREST.${identity}_${overlap}.tsv
rm all_SRG_related_ROIs.renamed.NR/EXACTMATCHES.${identity}_${overlap}.tsv
#index the insertion file in order to the get the contig sizes and names easily
samtools faidx all_SRG_related_ROIs.renamed.fa
#take the all insertions and blast each one seperately against the whole set
cat all_SRG_related_ROIs.renamed.fa.fai | while read line
do
contig=$( echo $line | awk '{print $1}' )
size=$( echo $line | awk '{print $2}' )
samtools faidx all_SRG_related_ROIs.renamed.fa "${contig}" > all_SRG_related_ROIs.renamed.NR/temp.fa
makeblastdb -in all_SRG_related_ROIs.renamed.NR/temp.fa -input_type fasta -parse_seqids -dbtype nucl
blastn -db all_SRG_related_ROIs.renamed.NR/temp.fa  -query all_SRG_related_ROIs.renamed.fa -outfmt 6 | awk -v contig="$contig" '{if($1 != contig) print}' >  all_SRG_related_ROIs.renamed.NR/${contig}.BLAST_raw.tsv
#take only alignments greater than the identity value and 5kb, rearrange inverted alignments, sort them and then generate non-redundant single alignments that are joint if there is a distance between them <= the gap value
cat all_SRG_related_ROIs.renamed.NR/${contig}.BLAST_raw.tsv | awk -v identity=".$identity" -v gap="$gap" '{if($3 > identity && $4 > gap) print $0 }' | awk '{if($9 < $10) {print $1"\t"$9"\t"$10} else {print $1"\t"$10"\t"$9}}' | sort -k1,1 -k2,2n | bedtools merge -d ${gap} > all_SRG_related_ROIs.renamed.NR/${contig}.BLAST_raw.merged2k.bed
#remove file incase the same output directory is being used again
rm all_SRG_related_ROIs.renamed.NR/${contig}.BLAST_raw.merged2k.proportion.bed
#get the size of the corresponding contig in order to calulate which contig is covered more by the alignment
cat all_SRG_related_ROIs.renamed.NR/${contig}.BLAST_raw.merged2k.bed | cut -f1 | sort -u | while read contig2
do
size2=$( cat all_SRG_related_ROIs.renamed.fa.fai | awk -v contig2="$contig2" '{if($1 == contig2) print $2}' )
cat all_SRG_related_ROIs.renamed.NR/${contig}.BLAST_raw.merged2k.bed | awk -v contig2="$contig2" -v size2="$size2" -v contig="$contig" -v size="$size" '{if($1 == contig2) print $2"\t"$3"\t"($3-$2)+1"\t"contig"\t"size"\t"(($3-$2)+1)/size"\t"contig2"\t"size2"\t"(($3-$2)+1)/size2}' >> all_SRG_related_ROIs.renamed.NR/${contig}.BLAST_raw.merged2k.proportion.bed
done
#take all alignments of interest, in that they are greater than the overlap value and that the contig in the 4th column is covered more by the alignment than the other contig
cat all_SRG_related_ROIs.renamed.NR/${contig}.BLAST_raw.merged2k.proportion.bed  | awk -v contig="$contig" -v size="$size" -v overlap=".$overlap" '{if($6 > overlap && $6 > $9)print $0}'  >> all_SRG_related_ROIs.renamed.NR/OFINTEREST.${identity}_${overlap}.tsv
cat all_SRG_related_ROIs.renamed.NR/${contig}.BLAST_raw.merged2k.proportion.bed  | awk -v contig="$contig" -v size="$size" -v overlap=".$overlap" '{if($6 > overlap && $6 == $9)print $0}'  >> all_SRG_related_ROIs.renamed.NR/EXACTMATCHES.${identity}_${overlap}.tsv
done
##have some exact matches and need to filter at least one of the matches
cat all_SRG_related_ROIs.renamed.NR/EXACTMATCHES.${identity}_${overlap}.tsv | awk '{print $4"\t"$7}' | while read line; do echo $line | awk '{print $1"\n"$2}' | sort | awk '{if(NR == 1){line=$1} else {line=line";"$1}} END {print line}' ; done | sort | uniq -c | awk '{if($1 == 2) print $2}' | awk -F ";" '{print ".\t.\t.\t"$1}'  >> all_SRG_related_ROIs.renamed.NR/OFINTEREST.${identity}_${overlap}.tsv
#get a list of the contigs to remove
list=$( cat all_SRG_related_ROIs.renamed.NR/OFINTEREST.${identity}_${overlap}.tsv | cut -f4 | sort -u )
#take the opposite for those to keep
listkeep=$(  cut -f1 all_SRG_related_ROIs.renamed.fa.fai | grep -v "${list}" )
#extract only those contigs to keep
samtools faidx all_SRG_related_ROIs.renamed.fa ${listkeep} > all_SRG_related_ROIs.renamed.NR.fa

##now we have out non redundant list of events we can extract the list of remaining HTRs and grab the SRGs related to them
samtools faidx all_SRG_related_ROIs.renamed.NR.fa
echo "nonredundant_name;original_name;annotation_name;size;SRGs;classification" | tr ';' '\t' > all_SRG_related_ROIs.renamed.NR.size_SRGclass.tsv
cat all_SRG_related_ROIs.renamed.NR.fa.fai | while read line
do
name=$( echo $line | awk '{print $1}' )
name2=$( cat all_SRG_related_ROIs.name_mod_association.tsv | awk -v name="$name" '{if($2 == name) print $1}' )
name3=$( cat **/ROI_annotation_${minsizeinsertion2}kb_min/name_mod_association.tsv | awk -v name="$name2" '{if($1 == name) print $2}' )
size=$( echo $line | awk '{print $2}' )
class=$( cat **/ROI_annotation_${minsizeinsertion2}kb_min/*.LR_combined.SRGs_positions.tsv | awk -v name="$name3" '{if($1 == name) genes=genes";"$6} END{print genes}' | sed "s/\t;/\t/g" )
class2=$( echo $class | awk '{if($1 ~ "MYB" && $1 ~ "DUF3435") {print "DUF3435_MYB/SANT"} else if($1 ~ "MYB" && $1 !~ "DUF3435") {print "MYB/SANT"} else if($1 !~ "MYB" && $1 ~ "DUF3435") {print "DUF3435"} else {print "Neither"}}' )
echo $name $name2 $name3 $size $class $class2 | tr ' ' '\t'
done >> all_SRG_related_ROIs.renamed.NR.size_SRGclass.tsv

###in order to show that DUF3435s and MYB/SANT genes sit at opposite ends to the other, we can calculate the relative distance of non-captain SRGs to DUF3435s
cat **/ROI_annotation_${minsizeinsertion2}kb_min/*.LR_combined.SRGs_positions.tsv | head -n1 | sed 's/contig/ROI/g' | awk '{print $0"\tcaptain_position\tabsolute_distance\trelative_distance"}' > all_SRG_related_ROIs.renamed.NR.noncaptains_distance_captains.tsv
cat all_SRG_related_ROIs.renamed.NR.size_SRGclass.tsv | awk '{if($6 ~ "DUF3435") print}' | while read line
do HTR=$( echo $line | awk '{print $3}' )
size=$( cat **/ROI_annotation_${minsizeinsertion2}kb_min/*.LR_combined.SRGs_positions.contigs_size.tsv | awk -v HTR="$HTR" '{if($1 == HTR) print $3}'  )
cat **/ROI_annotation_${minsizeinsertion2}kb_min/*.LR_combined.SRGs_positions.tsv | awk -v HTR="$HTR" '{if($1 == HTR && $6 == "DUF3435") print}' | while read captains
do
position=$( echo $captains | awk '{print ($3+$2)/2}' )
cat **/ROI_annotation_${minsizeinsertion2}kb_min/*.LR_combined.SRGs_positions.tsv | awk -v HTR="$HTR" '{if($1 == HTR && $6 != "DUF3435") print}' | while read noncaptains
do
echo $noncaptains" "$position | tr ' ' '\t' | awk '{print $0"\t"$7-(($3+$2)/2)}' | sed 's/-//g' | awk -v size="$size" '{print $0"\t"$NF/size}'
done
done
done >> all_SRG_related_ROIs.renamed.NR.noncaptains_distance_captains.tsv
