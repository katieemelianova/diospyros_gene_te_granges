
library(GenomicRanges)
library(gggenes)
library(plyranges)
source("/Users/katieemelianova/Desktop/Diospyros/diospyros_R_functions/diospyros_genome_annotation.R")



#########################################################
#                read in gene annotation                #
#########################################################

gene_granges <- get_braker_annotation()

#########################################################
#                read in EDTA annotation                #
#########################################################

te_intact_granges <- get_edta_intact_annotation()

impolita_te_granges <- te_intact_granges$impolita_edta_intact %>%
  separate(annotation, sep=";", c("ID", "Name", "Classification", "sequence_ontology")) %>%
  filter(str_detect(Name, "^Name")) %>%     # keep only entries where the TE is the parent TE
  separate(Name, sep="=", c("ph1", "Name")) %>%
  dplyr::select(-ph1)

impolita_te_granges$seqname <- paste0("Scaffolds_", impolita_te_granges$seqname) # change back names to match braker contigs names, EDTA edits them to make them shorter




#########################################################
#              read in TEsorter annotation              #
#########################################################

tesorter_granges <- get_tesorter_annotation()

impolita_tesorter <- tesorter_granges$impolita_tesorter %>%
  separate(detail, sep="#", c("Name", "ph1")) %>%
  separate(Name, sep="=", c("ph2", "Name")) %>%
  mutate(Name=str_remove(Name, "_INT")) %>%
  dplyr::select(-c(ph1, ph2))

# get gene and TE density using this method:
# https://www.biostars.org/p/169171/



#########################################################
#              read in TEsorter annotation              #
#########################################################

impolita_classified_te_intact_granges <- left_join(impolita_te_granges, impolita_tesorter, by="Name", relationship = "many-to-many") %>% dplyr::select(seqname, method, feature, start, end, Name, Classification, superclass, clade) %>% unique() %>% GRanges()



#########################################################
#              read in TEsorter annotation              #
#########################################################


impolita_seqlengths <- read.table("/Users/katieemelianova/Desktop/Diospyros/IGVdata/impolita/impolita.seqlengths")
impolita_seqlengths <- Seqinfo(seqnames=impolita_seqlengths$V1, seqlengths=impolita_seqlengths$V2)
impolita_gene_granges <- gene_granges$impolita %>% GRanges() %>% plyranges::filter(feature == "transcript")
impolita_tiles <- tileGenome(impolita_seqlengths, tilewidth=10000, cut.last.tile.in.chrom=T)
impolita_tiles$total_tes = countOverlaps(impolita_tiles, impolita_classified_te_intact_granges)
impolita_tiles$total_genes = countOverlaps(impolita_tiles, impolita_gene_granges)




impolita_clade_counts <- impolita_classified_te_intact_granges$clade %>% table() %>% data.frame() %>% set_colnames(c("clade", "freq"))
test_hidensity <- impolita_tiles %>% plyranges::filter(total_tes >= 1 & total_genes >= 2) %>% findOverlaps(impolita_classified_te_intact_granges)
test_lodensity <- impolita_tiles %>% plyranges::filter(total_tes >= 1 & total_genes == 0) %>% findOverlaps(impolita_classified_te_intact_granges)

# I think that some percentages are over 100% because of some TEs spanning multiple windows
# need to sanity check this

rbind(impolita_classified_te_intact_granges[test_hidensity@to %>% unique()] %>% data.frame() %>% pull(clade) %>% table() %>% data.frame() %>% mutate(density="high") %>% set_colnames(c("clade", "freq", "density")) %>% left_join(impolita_clade_counts, by="clade"),
      impolita_classified_te_intact_granges[test_lodensity@to %>% unique()] %>% data.frame() %>% pull(clade) %>% table() %>% data.frame() %>% mutate(density="low") %>% set_colnames(c("clade", "freq", "density")) %>% left_join(impolita_clade_counts, by="clade")) %>%
  mutate(percent_freq=(freq.x/freq.y) * 100) %>%
  filter(!(clade %in% c("Ty3_gypsy", "Ty1_copia", "chromo-unclass"))) %>%
  ggplot(aes(x=clade, y=percent_freq, fill=density)) +
  geom_bar(stat="identity")







summary(impolita_tiles$total_tes)
summary(impolita_tiles$total_genes)

# ask which TE classes are found how frequently in the same window as genes



mp<-readMappings("/Users/katieemelianova/Desktop/Diospyros/diospyros_gene_te_overlap/impolita_topGO_annotation.txt")
impolita_annotation <- read_delim("/Users/katieemelianova/Desktop/Diospyros/diospyros_gene_te_overlap/impolita_topGO_annotation.txt")

test <- findOverlaps(impolita_gene_granges, impolita_tiles %>% plyranges::filter(total_genes > 5 & total_tes > 5))
test2 <- impolita_gene_granges[test@from]$annotation
test2 %>% get_enriched_terms(mp, return_sample_GOData=TRUE)




impolita_annotation %>% filter(ID %in% test2) %>% pull(GO)


imp_overlaps <- findOverlaps(gene_granges$impolita, te_intact_granges$impolita)
rev_overlaps <- findOverlaps(gene_granges$revolutissima, te_intact_granges$revolutissima)
vie_overlaps <- findOverlaps(gene_granges$vieillardii, te_intact_granges$vieillardii)

feature_overlap_species <- rbind(gene_granges$impolita[imp_overlaps@from]$feature %>% table() %>% data.frame() %>% filter(!(. == "mRNA")) %>% mutate(species = "D. impolita"),
gene_granges$revolutissima[rev_overlaps@from]$feature %>% table() %>% data.frame() %>% filter(!(. == "mRNA")) %>% mutate(species = "D. revolutissima"),
gene_granges$vieillardii[vie_overlaps@from]$feature %>% table() %>% data.frame() %>% filter(!(. == "mRNA")) %>% mutate(species = "D. vieillardii")) %>%
  set_colnames(c("feature", "freq", "species")) %>%
  filter(feature %in% c("CDS", "exon", "intron", "start_codon", "stop_codon") & !(species == "D. yahouensis"))

feature_overlap_species$species <- factor(feature_overlap_species$species, levels=c("D. sandwicensis", "D. vieillardii", "D. pancheri", "D. yahouensis", "D. revolutissima", "D. impolita"))

ggplot(data=feature_overlap_species, aes(x=feature, y=(log(freq)), fill=species)) +
  geom_bar(stat="identity", position = "dodge")






















