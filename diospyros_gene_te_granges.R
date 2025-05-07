
library(GenomicRanges)
library(gggenes)
library(plyranges)




#########################################################
#                read in gene annotation                #
#########################################################

read_genes <- function(gff){
  genes_gff <- read_delim(gff, col_names = FALSE) %>%
    set_colnames(c("seqname", "method", "feature", "start", "end", "ph1", "ph2", "ph3", "annotation")) %>% 
    GRanges()
  return(genes_gff)
}

gene_files <- c("/Users/katieemelianova/Desktop/Diospyros/IGVdata/impolita/impolita_braker.gtf",
                "/Users/katieemelianova/Desktop/Diospyros/IGVdata/sandwicensis/sandwicensis_braker.gtf",
                "/Users/katieemelianova/Desktop/Diospyros/IGVdata/pancheri/pancheri_braker.gtf",
                "/Users/katieemelianova/Desktop/Diospyros/IGVdata/vieillardii/vieillardii_braker.gtf",
                "/Users/katieemelianova/Desktop/Diospyros/IGVdata/revolutissima/revolutissima_braker.gtf",
                "/Users/katieemelianova/Desktop/Diospyros/IGVdata/yahouensis/yahouensis_braker.gtf")


gene_granges <- gene_files %>%
  lapply(read_genes) %>% 
  set_names(c("impolita", "sandwicensis", "pancheri", "vieillardii", "revolutissima", "yahouensis"))




#########################################################
#                read in gene annotation                #
#########################################################

read_tes <- function(gff){
  tes_gff <- read_delim(gff) %>%
    set_colnames(c("seqname", "method", "feature", "start", "end", "ph1", "ph2", "ph3", "annotation")) %>%
    filter(feature == "repeat_region") %>%
    GRanges()
  return(tes_gff)
}

te_intact_files <- c("/Users/katieemelianova/Desktop/Diospyros/IGVdata/impolita/impolita.fasta.mod.EDTA.intact.noheader.seqnames.gff3",
                     "/Users/katieemelianova/Desktop/Diospyros/IGVdata/pancheri/pancheri.fasta.mod.EDTA.intact.noheader.gff3",
                     "/Users/katieemelianova/Desktop/Diospyros/IGVdata/revolutissima/revolutissima.fasta.mod.EDTA.intact.noheader.gff3",
                     "/Users/katieemelianova/Desktop/Diospyros/IGVdata/sandwicensis/sandwicensis.fasta.mod.EDTA.intact.noheader.gff3",
                     "/Users/katieemelianova/Desktop/Diospyros/IGVdata/vieillardii/vieillardii.fasta.mod.EDTA.intact.noheader.gff3",
                     "/Users/katieemelianova/Desktop/Diospyros/IGVdata/yahouensis/yahouensis.fasta.mod.EDTA.intact.noheader.gff3")

te_intact_granges <- te_intact_files %>%
  lapply(read_tes) %>% 
  set_names(c("impolita", "pancheri", "revolutissima", "sandwicensis", "vieillardii", "yahouensis"))

# get gene and TE density using this method:
# https://www.biostars.org/p/169171/

impolita_seqlengths <- read.table("/Users/katieemelianova/Desktop/Diospyros/IGVdata/impolita/impolita.seqlengths")
impolita_seqlengths <- Seqinfo(seqnames=impolita_seqlengths$V1, seqlengths=impolita_seqlengths$V2)
impolita_te_granges <- te_intact_granges$impolita
impolita_gene_granges <- gene_granges$impolita %>% plyranges::filter(feature == "transcript")
impolita_tiles <- tileGenome(impolita_seqlengths, tilewidth=100000, cut.last.tile.in.chrom=T)
impolita_tiles$total_tes = countOverlaps(impolita_tiles, impolita_te_granges)
impolita_tiles$total_genes = countOverlaps(impolita_tiles, impolita_gene_granges)

plot(impolita_tiles$total_tes, impolita_tiles$total_genes)


impolita_tiles$total_tes %>% summary()
impolita_tiles$total_genes %>% summary()



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






















