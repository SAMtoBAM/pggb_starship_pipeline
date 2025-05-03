###this is a modifiable Rscript for generating Starship alignments
##the basic skeleton will remain the same and the same features can be swapped in for each variable

##only need gggenomes
library(gggenomes)
library(ggnewscale)

##need four features

##First feature
##first a bed file with just the length and size of the full alignment regions
bed=read.csv("~/projects/Penicillium/genome_graphs/SLR_plots/between_SLRs_plus/cluster_IBT27055:HTR44.SRGs.50kb_flank.bed", sep='\t', header=T)
##just modify some headers for downstream handling
##make another header, copy of contig, but called seq_id
bed$seq_id = bed$contig
##also length using the end position
bed$length = bed$end
#bed$length = bed$end - bed$start

##second feature
##a bed file of just the Starship-like region coordinates i.e. the above bed file without the flanking regions
SLRbed=read.csv("~/projects/Penicillium/genome_graphs/SLR_plots/between_SLRs_plus/cluster_IBT27055:HTR44.SRGs.50kb_flank.SRG.bed", sep='\t', header=T)
SLRbed$seq_id = SLRbed$contig
SLRbed$length = SLRbed$end-SLRbed$start

##third feature
##a bed file with the genes annotated within the Starship-like regions (coordinates modified due to the flanking regions added)
genes=read.csv("~/projects/Penicillium/genome_graphs/SLR_plots/between_SLRs_plus/cluster_IBT27055:HTR44.SRGs.50kb_flank.genes.bed", sep='\t', header=T)
genes$seq_id = genes$contig
genes$length= genes$end-genes$start
genes$strand=genes$sense

##fourth feature
##the nucmer all-v-all alignment converted to paf format
links=read_links("~/projects/Penicillium/genome_graphs/SLR_plots/between_SLRs_plus/cluster_IBT27055:HTR44.SRGs.50kb_flank.nucmer.paf")




##the actual plot
## only selecting alignments greater than 5kb and with greater then _(% identity)
gggenomes(genes=genes, seqs=bed, feat=SLRbed, links=subset(links, map_length > 10000 & map_match/map_length > 0.8 & seq_id != seq_id2), adjacent_only = F) %>%
  pick() %>%
  flip() %>%
  gggenomes::sync() +
  geom_link(aes(fill=map_match/map_length) ,colour="black", alpha=0.5, offset = 0.05, size=0.1 )+
  scale_fill_gradientn(colours=c("grey","orange", "orange4"), name ="Identity", labels=c(.80,0.90,1), breaks=c(0.80,0.90,1), limits = c(0.8, 1))+
  new_scale_fill()+
  geom_seq(linewidth = 0.5)+
  geom_feat(color="lightblue", alpha=.6, linewidth=3)+
  geom_gene(aes(fill=label), stroke=0.1, colour="black", shape = 3)+
  geom_seq_label()+
  geom_seq_label(aes(label=fullname), nudge_y = -.25)+
  geom_gene_tag(aes(label=label), size = 2, nudge_y=0.1, check_overlap = FALSE)+
  scale_fill_manual(values = c("red","blue","purple","yellow","orange","darkgreen"), breaks=c("DUF3435","MYB/SANT", "DUF3723","PLP","CRG", "DUF6066"))+
  theme(legend.position="bottom", legend.box = "horizontal")

##save plot as variable so can save it
plot=gggenomes(genes=genes, seqs=bed, feat=SLRbed, links=subset(links, map_length > 8000 & map_match/map_length > 0.85 & seq_id != seq_id2), adjacent_only = F) %>%
  pick() %>%
  flip() %>%
  gggenomes::sync() +
  geom_link(aes(fill=map_match/map_length) ,colour="black", alpha=0.5, offset = 0.05, size=0.1 )+
  scale_fill_gradientn(colours=c("grey","orange", "orange4"), name ="Identity", labels=c(.80,0.90,1), breaks=c(0.80,0.90,1), limits = c(0.8, 1))+
  new_scale_fill()+
  geom_seq(linewidth = 0.5)+
  geom_feat(color="lightblue", alpha=.6, linewidth=3)+
  geom_gene(aes(fill=label), stroke=0.1, colour="black", shape = 3)+
  geom_seq_label()+
  geom_seq_label(aes(label=fullname), nudge_y = -.25)+
  geom_gene_tag(aes(label=label), size = 2, nudge_y=0.1, check_overlap = FALSE)+
  scale_fill_manual(values = c("red","blue","purple","yellow","orange","darkgreen"), breaks=c("DUF3435","MYB/SANT", "DUF3723","PLP","CRG", "DUF6066"))+
  theme(legend.position="bottom", legend.box = "horizontal")


ggsave(filename="~/projects/Penicillium/genome_graphs/SLR_plots/between_SLRs_plus/cluster_IBT27055:HTR44.SRGs.50kb_flank.pdf", 
       plot = plot, 
       device = cairo_pdf, 
       width = 210, 
       height = 297, 
       units = "mm")

plot2=gggenomes(genes=genes, seqs=bed, feat=SLRbed, links=subset(links, map_match/map_length > 0.8 & seq_id != seq_id2)) %>%
  pick() %>%
  flip() %>%
  gggenomes::sync() +
  geom_link(aes(fill=map_match/map_length) ,colour="black", alpha=0.5, offset = 0.05, size=0.1 )+
  scale_fill_gradientn(colours=c("grey","orange", "orange4"), name ="Identity", labels=c(.80,0.90,1), breaks=c(0.80,0.90,1), limits = c(0.8, 1))+
  new_scale_fill()+
  geom_seq(linewidth = 0.5)+
  geom_feat(color="lightblue", alpha=.6, linewidth=3)+
  geom_gene(aes(fill=label), stroke=0.1, colour="black", shape = 3)+
  geom_seq_label()+
  geom_seq_label(aes(label=fullname), nudge_y = -.25)+
  geom_gene_tag(aes(label=label), size = 2, nudge_y=0.1, check_overlap = FALSE)+
  scale_fill_manual(values = c("red","blue","purple","yellow","orange","darkgreen"), breaks=c("DUF3435","MYB/SANT", "DUF3723","PLP","CRG", "DUF6066"))+
  theme(legend.position="bottom", legend.box = "horizontal")


ggsave(filename="~/projects/Penicillium/genome_graphs/SLR_plots/between_SLRs_plus/cluster_IBT27055:HTR44.SRGs.50kb_flank.adjacent_only.all_alignments.pdf", 
       plot = plot2, 
       device = cairo_pdf, 
       width = 210, 
       height = 297, 
       units = "mm")


plot3=gggenomes(genes=genes, seqs=bed, feat=SLRbed, links=subset(links, map_length > 8000 & map_match/map_length > 0.8 & seq_id != seq_id2)) %>%
  pick() %>%
  flip() %>%
  gggenomes::sync() +
  geom_link(aes(fill=map_match/map_length) ,colour="black", alpha=0.5, offset = 0.05, size=0.1 )+
  scale_fill_gradientn(colours=c("grey","orange", "orange4"), name ="Identity", labels=c(.80,0.90,1), breaks=c(0.80,0.90,1), limits = c(0.8, 1))+
  new_scale_fill()+
  geom_seq(linewidth = 0.5)+
  geom_feat(color="lightblue", alpha=.6, linewidth=3)+
  geom_gene(aes(fill=label), stroke=0.1, colour="black", shape = 3)+
  geom_seq_label()+
  geom_seq_label(aes(label=fullname), nudge_y = -.25)+
  geom_gene_tag(aes(label=label), size = 2, nudge_y=0.1, check_overlap = FALSE)+
  scale_fill_manual(values = c("red","blue","purple","yellow","orange","darkgreen"), breaks=c("DUF3435","MYB/SANT", "DUF3723","PLP","CRG", "DUF6066"))+
  theme(legend.position="bottom", legend.box = "horizontal")


ggsave(filename="~/projects/Penicillium/genome_graphs/SLR_plots/between_SLRs_plus/cluster_IBT27055:HTR44.SRGs.50kb_flank.adjacent_only.pdf", 
       plot = plot3, 
       device = cairo_pdf, 
       width = 210, 
       height = 297, 
       units = "mm")
