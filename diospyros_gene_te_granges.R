
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

# change back names to match braker contigs names, EDTA edits them to make them shorter (only an issue in impolita)
impolita_te_granges$seqname <- paste0("Scaffolds_", impolita_te_granges$seqname) 

revolutissima_te_granges <- te_intact_granges$revolutissima_edta_intact %>%
  separate(annotation, sep=";", c("ID", "Name", "Classification", "sequence_ontology")) %>%
  filter(str_detect(Name, "^Name")) %>%     # keep only entries where the TE is the parent TE
  separate(Name, sep="=", c("ph1", "Name")) %>%
  dplyr::select(-ph1)

vieillardii_te_granges <- te_intact_granges$vieillardii_edta_intact %>%
  separate(annotation, sep=";", c("ID", "Name", "Classification", "sequence_ontology")) %>%
  filter(str_detect(Name, "^Name")) %>%     # keep only entries where the TE is the parent TE
  separate(Name, sep="=", c("ph1", "Name")) %>%
  dplyr::select(-ph1)

#########################################################
#              read in TEsorter annotation              #
#########################################################

tesorter_granges <- get_tesorter_annotation()

impolita_tesorter <- tesorter_granges$impolita_tesorter %>%
  separate(detail, sep="#", c("Name", "ph1")) %>%
  separate(Name, sep="=", c("ph2", "Name")) %>%
  mutate(Name=str_remove(Name, "_INT")) %>%
  dplyr::select(-c(ph1, ph2))

revolutissima_tesorter <- tesorter_granges$revolutissima_tesorter %>%
  separate(detail, sep="#", c("Name", "ph1")) %>%
  separate(Name, sep="=", c("ph2", "Name")) %>%
  mutate(Name=str_remove(Name, "_INT")) %>%
  dplyr::select(-c(ph1, ph2))

vieillardii_tesorter <- tesorter_granges$vieillardii_tesorter %>%
  separate(detail, sep="#", c("Name", "ph1")) %>%
  separate(Name, sep="=", c("ph2", "Name")) %>%
  mutate(Name=str_remove(Name, "_INT")) %>%
  dplyr::select(-c(ph1, ph2))
# get gene and TE density using this method:
# https://www.biostars.org/p/169171/



#########################################################
#           join EDTA and TEsorter annotation           #
#########################################################

impolita_classified_te_intact_granges <- left_join(impolita_te_granges, impolita_tesorter, by="Name", relationship = "many-to-many") %>% dplyr::select(seqname, method, feature, start, end, Name, Classification, superclass, clade) %>% unique() %>% GRanges()
revolutissima_classified_te_intact_granges <- left_join(revolutissima_te_granges, revolutissima_tesorter, by="Name", relationship = "many-to-many") %>% dplyr::select(seqname, method, feature, start, end, Name, Classification, superclass, clade) %>% unique() %>% GRanges()
vieillardii_classified_te_intact_granges <- left_join(vieillardii_te_granges, vieillardii_tesorter, by="Name", relationship = "many-to-many") %>% dplyr::select(seqname, method, feature, start, end, Name, Classification, superclass, clade) %>% unique() %>% GRanges()



#########################################################
#      tile genomes and plot gene and TE density        #
#########################################################

tile_genome <- function(seqlengths, gene_grange, classified_te_grange){
  # tile genome using seqlengths generated with bioawk
  seqlengths <- read.table(seqlengths)
  seqlengths <- Seqinfo(seqnames=seqlengths$V1, seqlengths=seqlengths$V2)
  gene_grange <- gene_grange %>% GRanges() %>% plyranges::filter(feature == "transcript")
  tiles <- tileGenome(seqlengths, tilewidth=1000000, cut.last.tile.in.chrom=T)
  # count overlaps of tiles with genes and TEs
  tiles$total_tes = countOverlaps(tiles, classified_te_grange, ignore.strand=TRUE)
  tiles$total_genes = countOverlaps(tiles, gene_grange, ignore.strand=TRUE)
  # get clade counts- these are used to calculate the percent of each TE clade in high and low gene density regions
  clade_counts <- classified_te_grange$clade %>% table() %>% data.frame() %>% set_colnames(c("clade", "freq"))
  # get high and low density tiles, get the overlap with the annotated TE object and get the clades present in high and low density regions
  hidensity <- tiles %>% plyranges::filter(total_tes >= 1 & total_genes >= 20) %>% findOverlaps(classified_te_grange)
  lodensity <- tiles %>% plyranges::filter(total_tes >= 1 & total_genes == 0) %>% findOverlaps(classified_te_grange)
  totdensity <- clade_counts %>% mutate(density="total", freq.y=freq) %>% set_colnames(c("clade", "freq.x", "density", "freq.y"))
  to_return <- rbind(classified_te_grange[hidensity@to %>% unique()] %>% data.frame() %>% pull(clade) %>% table() %>% data.frame() %>% mutate(density="high") %>% set_colnames(c("clade", "freq", "density")) %>% left_join(clade_counts, by="clade"),
                     classified_te_grange[lodensity@to %>% unique()] %>% data.frame() %>% pull(clade) %>% table() %>% data.frame() %>% mutate(density="low") %>% set_colnames(c("clade", "freq", "density")) %>% left_join(clade_counts, by="clade")) %>%
    rbind(totdensity) %>%
    mutate(percent_freq=(freq.x/freq.y) * 100) %>%
    filter(!(clade %in% c("Ty3_gypsy", "Ty1_copia", "chromo-unclass")))
  
  return(to_return)
}

impolita_tiles <- tile_genome("/Users/katieemelianova/Desktop/Diospyros/IGVdata/impolita/impolita.seqlengths",
            gene_granges$impolita_braker,
            impolita_classified_te_intact_granges) %>%
  mutate(species="impolita")

revolutissima_tiles <- tile_genome("/Users/katieemelianova/Desktop/Diospyros/IGVdata/revolutissima/revolutissima.seqlengths",
                              gene_granges$revolutissima_braker,
                              revolutissima_classified_te_intact_granges) %>%
  mutate(species="revolutissima")

vieillardii_tiles <- tile_genome("/Users/katieemelianova/Desktop/Diospyros/IGVdata/vieillardii/vieillardii.seqlengths",
                                 gene_granges$vieillardii_braker,
                                 vieillardii_classified_te_intact_granges) %>%
  mutate(species="vieillardii")


# find clades which are not present in all three species, we wont plot these
include_clades <- rbind(impolita_tiles,
      revolutissima_tiles,
      vieillardii_tiles)%>%
  dplyr::select(-c(freq.x, freq.y)) %>%
  melt() %>%
  dplyr::count(clade) %>%
  tibble() %>%
  filter(n >=3) %>%
  pull(clade) %>%
  as.character()





rbind(impolita_tiles,
        revolutissima_tiles,
        vieillardii_tiles) %>%
  mutate(density=factor(density, levels=c("low", "high", "total"))) %>% 
  arrange(desc(density)) %>%
  filter(clade %in% include_clades) %>%
  ggplot(aes(x=clade, y=sqrt(freq.x), fill=density)) +
  geom_bar(stat="identity", position = "identity", alpha=0.7) +
  facet_wrap(~species, ncol = 1) +
  scale_fill_manual(values=c("blue", "red", "grey67"))

# have a look to see how many TEs span more than one window
# I increased window size to 1MB because 
findOverlaps(revolutissima_tiles, revolutissima_classified_te_intact_granges) %>% 
  data.frame() %>% 
  pull(subjectHits) %>%
  table() %>% 
  sort(decreasing = TRUE) %>% 
  data.frame()










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






















