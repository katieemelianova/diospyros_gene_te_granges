
library(GenomicRanges)
library(gggenes)




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




###########################################################################
#                        Plot Annexin OG0000336                           #
###########################################################################

# vieillardii "g1861"
# revolutissima "g11314"
# impolita "g4621"


revolutissima_g11314_gene <- gene_granges$revolutissima[gene_granges$revolutissima$annotation == "g11314"]
vieillardii_g1861_gene <- gene_granges$vieillardii[gene_granges$vieillardii$annotation == "g1861"]
impolita_g4621_gene <- gene_granges$impolita[gene_granges$impolita$annotation == "g4621"]

revolutissima_g11314_intact_TE <- findOverlaps(revolutissima_g11314_gene, te_intact_granges$revolutissima,  maxgap = 35000)
vieillardii_g1861_intact_TE <- findOverlaps(vieillardii_g1861_gene, te_intact_granges$vieillardii,  maxgap = 35000)
impolita_g4621_intact_TE <- findOverlaps(impolita_g4621_gene, te_intact_granges$impolita,  maxgap = 35000)

revolutissima_g11314_intact_TE <- te_intact_granges$revolutissima[revolutissima_g11314_intact_TE@to]
vieillardii_g1861_intact_TE <- te_intact_granges$vieillardii[vieillardii_g1861_intact_TE@to]
impolita_g4621_intact_TE <- te_intact_granges$impolita[impolita_g4621_intact_TE@to]

test_te <- rbind(revolutissima_g11314_intact_TE %>% data.frame() %>% mutate(gene="TE"), 
      vieillardii_g1861_intact_TE %>% data.frame() %>% mutate(gene="TE"), 
      impolita_g4621_intact_TE %>% data.frame()%>% mutate(gene="TE")) %>% 
  dplyr::select(c(seqnames, gene, start, end, ph2))


test_gene <- rbind(revolutissima_g11314_gene %>% data.frame() %>% mutate(gene="Annexin"),
vieillardii_g1861_gene %>% data.frame() %>% mutate(gene="Annexin"),
impolita_g4621_gene %>% data.frame() %>% mutate(gene="Annexin")) %>% 
  dplyr::select(c(seqnames, gene, start, end, ph2))

test <- rbind(test_te, test_gene) %>% mutate(species = case_when(seqnames == "ptg000017l" ~ "revolutissima",
                                                                 seqnames == "ptg000002l" ~ "vieillardii",
                                                                 seqnames == "Scaffolds_1151" ~ "impolita"))

dummies <- make_alignment_dummies(
  test,
  aes(xmin = start, xmax = end, y = species, id = gene),
  on = "Annexin"
)

pdf("annexin_OG0000336_geneplot.pdf", height=5.5, width=9)
ggplot2::ggplot(test, ggplot2::aes(xmin = start, xmax = end,
                                            y = species, fill = gene, label = gene)) +
  geom_gene_arrow(arrow_body_height = grid::unit(10, "mm"),
                  arrowhead_height = grid::unit(12, "mm")) +
  #geom_gene_label(height = grid::unit(6, "mm"), grow = TRUE) +
  ggplot2::facet_wrap(~ species, ncol = 1, scales = "free") +
  theme_genes() +
  theme(legend.text = element_text(size=20),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=14),
        axis.title = element_text(size=20),
        legend.title = element_blank()) +
  geom_blank(data = dummies) +
  ylab("Species")
dev.off()



###########################################################################
#                            Plot OG0000500                               #
###########################################################################


"vieillardii_braker_aa_g10152.t1"
"impolita_braker_aa_g22665.t1"
"revolutissima_braker_aa_g482.t1"

revolutissima_g482_gene <- gene_granges$revolutissima[gene_granges$revolutissima$annotation == "g482"]
vieillardii_g10152_gene <- gene_granges$vieillardii[gene_granges$vieillardii$annotation == "g10152"]
impolita_g22665_gene <- gene_granges$impolita[gene_granges$impolita$annotation == "g22665"]

revolutissima_g482_intact_TE <- findOverlaps(revolutissima_g482_gene, te_intact_granges$revolutissima,  maxgap = 50000)
vieillardii_g10152_intact_TE <- findOverlaps(vieillardii_g10152_gene, te_intact_granges$vieillardii,  maxgap = 50000)
impolita_g22665_intact_TE <- findOverlaps(impolita_g22665_gene, te_intact_granges$impolita,  maxgap = 50000)

revolutissima_g482_intact_TE <- te_intact_granges$revolutissima[revolutissima_g482_intact_TE@to]
vieillardii_g10152_intact_TE <- te_intact_granges$vieillardii[vieillardii_g10152_intact_TE@to]
impolita_g22665_intact_TE <- te_intact_granges$impolita[impolita_g22665_intact_TE@to]

OG0000500_te <- rbind(revolutissima_g482_intact_TE %>% data.frame() %>% mutate(gene="TE"), 
                      vieillardii_g10152_intact_TE %>% data.frame() %>% mutate(gene="TE"), 
                      impolita_g22665_intact_TE %>% data.frame()%>% mutate(gene="TE")) %>% 
  dplyr::select(c(seqnames, gene, start, end, ph2))


OG0000500_gene <- rbind(revolutissima_g482_gene %>% data.frame() %>% mutate(gene="PR4A"),
                        vieillardii_g10152_gene %>% data.frame() %>% mutate(gene="PR4A"),
                        impolita_g22665_gene %>% data.frame() %>% mutate(gene="PR4A")) %>% 
  dplyr::select(c(seqnames, gene, start, end, ph2))

test_pr4a <- rbind(OG0000500_te, OG0000500_gene) %>% mutate(species = case_when(seqnames == "ptg000001l" ~ "revolutissima",
                                                                 seqnames == "ptg000010l" ~ "vieillardii",
                                                                 seqnames == "Scaffolds_739" ~ "impolita"))



dummies <- make_alignment_dummies(
  test_pr4a,
  aes(xmin = start, xmax = end, y = species, id = gene),
  on = "PR4A"
)

pdf("PR4A_OG0000500_geneplot.pdf", height=5.5, width=9)
ggplot2::ggplot(test_pr4a, ggplot2::aes(xmin = start, xmax = end,
                                   y = species, fill = gene, label = gene)) +
  geom_gene_arrow(arrow_body_height = grid::unit(10, "mm"),
                  arrowhead_height = grid::unit(12, "mm")) +
  #geom_gene_label(height = grid::unit(6, "mm"), grow = TRUE) +
  ggplot2::facet_wrap(~ species, ncol = 1, scales = "free") +
  geom_blank(data = dummies) +
  theme_genes() +
  theme(legend.text = element_text(size=20),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=14),
        axis.title = element_text(size=20),
        legend.title = element_blank()) +
  ylab("Species")
dev.off()



###########################################################################
#                            Plot OG0006346                               #
###########################################################################


"vieillardii_braker_aa_g594.t1"
"revolutissima_braker_aa_g14832.t1"
"impolita_braker_aa_g3611.t1"


revolutissima_gene <- gene_granges$revolutissima[gene_granges$revolutissima$annotation == "g14832"]
vieillardii_gene <- gene_granges$vieillardii[gene_granges$vieillardii$annotation == "g594"]
impolita_gene <- gene_granges$impolita[gene_granges$impolita$annotation == "g3611"]
pancheri_gene <- gene_granges$pancheri[gene_granges$pancheri$annotation == "g25218"]
yahouensis_gene <- gene_granges$yahouensis[gene_granges$yahouensis$annotation == "g9561"]


revolutissima_intact_TE <- findOverlaps(revolutissima_gene, te_intact_granges$revolutissima,  maxgap = 170000)
vieillardii_intact_TE <- findOverlaps(vieillardii_gene, te_intact_granges$vieillardii,  maxgap = 170000)
impolita_intact_TE <- findOverlaps(impolita_gene, te_intact_granges$impolita,  maxgap = 170000)
pancheri_intact_TE <- findOverlaps(pancheri_gene, te_intact_granges$pancheri,  maxgap = 170000)
yahouensis_intact_TE <- findOverlaps(yahouensis_gene, te_intact_granges$yahouensis,  maxgap = 170000)

revolutissima_intact_TE <- te_intact_granges$revolutissima[revolutissima_intact_TE@to]
vieillardii_intact_TE <- te_intact_granges$vieillardii[vieillardii_intact_TE@to]
impolita_intact_TE <- te_intact_granges$impolita[impolita_intact_TE@to]
pancheri_intact_TE <- te_intact_granges$pancheri[pancheri_intact_TE@to]
yahouensis_intact_TE <- te_intact_granges$yahouensis[yahouensis_intact_TE@to]


#OG0000500_te <- rbind(revolutissima_intact_TE %>% data.frame() %>% mutate(gene="TE"), 
#                      vieillardii_intact_TE %>% data.frame() %>% mutate(gene="TE"), 
#                      impolita_intact_TE %>% data.frame()%>% mutate(gene="TE"),
#                      pancheri_intact_TE %>% data.frame() %>% mutate(gene="TE"),
#                      yahouensis_intact_TE %>% data.frame() %>% mutate(gene="TE")) %>% 
#  dplyr::select(c(seqnames, gene, start, end, ph2))
#
#
#OG0000500_gene <- rbind(revolutissima_gene %>% data.frame() %>% mutate(gene="SRC2"),
#                        vieillardii_gene %>% data.frame() %>% mutate(gene="SRC2"),
#                        impolita_gene %>% data.frame() %>% mutate(gene="SRC2"),
#                        pancheri_gene %>% data.frame() %>% mutate(gene="SRC2"),
#                        yahouensis_gene %>% data.frame() %>% mutate(gene="SRC2")) %>% 
#  dplyr::select(c(seqnames, gene, start, end, ph2))
#
#test_scr2 <- rbind(OG0000500_te, OG0000500_gene) %>% mutate(species = case_when(seqnames == "ptg000024l" ~ "revolutissima",
#                                                                                seqnames == "ptg000001l" ~ "vieillardii",
#                                                                                seqnames == "Scaffolds_1130" ~ "impolita",
#                                                                                seqnames == "ptg000121l" ~ "pancheri",
#                                                                                seqnames == "yahouensis_tig00007032" ~ "yahouensis"))
#


OG0000500_te <- rbind(revolutissima_intact_TE %>% data.frame() %>% mutate(gene="TE"), 
                      vieillardii_intact_TE %>% data.frame() %>% mutate(gene="TE"), 
                      impolita_intact_TE %>% data.frame()%>% mutate(gene="TE")) %>% 
  dplyr::select(c(seqnames, gene, start, end, ph2))


OG0000500_gene <- rbind(revolutissima_gene %>% data.frame() %>% mutate(gene="SRC2"),
                        vieillardii_gene %>% data.frame() %>% mutate(gene="SRC2"),
                        impolita_gene %>% data.frame() %>% mutate(gene="SRC2")) %>% 
  dplyr::select(c(seqnames, gene, start, end, ph2))

test_scr2 <- rbind(OG0000500_te, OG0000500_gene) %>% mutate(species = case_when(seqnames == "ptg000024l" ~ "revolutissima",
                                                                                seqnames == "ptg000001l" ~ "vieillardii",
                                                                                seqnames == "Scaffolds_1130" ~ "impolita"))

dummies <- make_alignment_dummies(
  test_scr2,
  aes(xmin = start, xmax = end, y = species, id = gene),
  on = "SRC2"
)

pdf("SRC2_OG0006346_geneplot.pdf", height=5.5, width=9)
ggplot2::ggplot(test_scr2, ggplot2::aes(xmin = start, xmax = end,
                                        y = species, fill = gene, label = gene)) +
  geom_gene_arrow(arrow_body_height = grid::unit(10, "mm"),
                  arrowhead_height = grid::unit(12, "mm")) +
  #geom_gene_label(height = grid::unit(6, "mm"), grow = TRUE) +
  ggplot2::facet_wrap(~ species, ncol = 1, scales = "free") +
  geom_blank(data = dummies) +
  theme_genes() +
  theme(legend.text = element_text(size=20),
        axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=14),
        axis.title = element_text(size=20),
        legend.title = element_blank()) +
  ylab("Species")
dev.off()








imp_overlaps <- findOverlaps(gene_granges$impolita, te_intact_granges$impolita)
rev_overlaps <- findOverlaps(gene_granges$revolutissima, te_intact_granges$revolutissima)
san_overlaps <- findOverlaps(gene_granges$sandwicensis, te_intact_granges$sandwicensis)
pan_overlaps <- findOverlaps(gene_granges$pancheri, te_intact_granges$pancheri)
vie_overlaps <- findOverlaps(gene_granges$vieillardii, te_intact_granges$vieillardii)
yah_overlaps <- findOverlaps(gene_granges$yahouensis, te_intact_granges$yahouensis)

feature_overlap_species <- rbind(gene_granges$impolita[imp_overlaps@from]$feature %>% table() %>% data.frame() %>% filter(!(. == "mRNA")) %>% mutate(species = "D. impolita"),
gene_granges$revolutissima[rev_overlaps@from]$feature %>% table() %>% data.frame() %>% filter(!(. == "mRNA")) %>% mutate(species = "D. revolutissima"),
gene_granges$sandwicensis[san_overlaps@from]$feature %>% table() %>% data.frame() %>% filter(!(. == "mRNA")) %>% mutate(species = "D. sandwicensis"),
gene_granges$pancheri[pan_overlaps@from]$feature %>% table() %>% data.frame() %>% filter(!(. == "mRNA")) %>% mutate(species = "D. pancheri"),
gene_granges$vieillardii[vie_overlaps@from]$feature %>% table() %>% data.frame() %>% filter(!(. == "mRNA")) %>% mutate(species = "D. vieillardii"),
gene_granges$yahouensis[yah_overlaps@from]$feature %>% table() %>% data.frame() %>% filter(!(. == "mRNA")) %>% mutate(species = "D. yahouensis")) %>%
  set_colnames(c("feature", "freq", "species")) %>%
  filter(feature %in% c("CDS", "exon", "intron", "start_codon", "stop_codon") & !(species == "D. yahouensis"))

feature_overlap_species$species <- factor(feature_overlap_species$species, levels=c("D. sandwicensis", "D. vieillardii", "D. pancheri", "D. yahouensis", "D. revolutissima", "D. impolita"))

ggplot(data=feature_overlap_species, aes(x=feature, y=(freq), fill=species)) +
  geom_bar(stat="identity", position = "dodge")



rev_overlaps

###







#########################################################
#                 get promoter regions                  #
#########################################################

impolita_promoters <- promoters(gene_granges$impolita)
revolutissima_promoters <- promoters(gene_granges$revolutissima)




impolita_te_promoter <- impolita_promoters[overlapsAny(impolita_promoters, te_intact_granges$impolita)]
revolutissima_te_promoter <- revolutissima_promoters[overlapsAny(revolutissima_promoters, te_intact_granges$revolutissima)]

imp_rev_de_orthogroups <- orthogroups_long %>% 
  filter(species == "D. vieillardii" & gene %in% (rev_imp$results %>% data.frame() %>% 
                                                    filter(log2FoldChange > 1.5 & padj < 0.05) %>% 
                                                    rownames() %>% str_split_i("\\.", 1))) %>% pull(Orthogroup)

impolita_de_orthologs <- orthogroups_long %>% filter(Orthogroup %in% imp_rev_de_orthogroups & species == "D. impolita") %>% drop_na()
revolutissima_de_orthologs <- orthogroups_long %>% filter(Orthogroup %in% imp_rev_de_orthogroups & species == "D. revolutissima") %>% drop_na()


imp_og <- impolita_de_orthologs %>% filter(gene %in% intersect(impolita_de_orthologs$gene, impolita_te_promoter$annotation)) %>% pull(Orthogroup)
rev_og<- revolutissima_de_orthologs %>% filter(gene %in% intersect(revolutissima_de_orthologs$gene, revolutissima_te_promoter$annotation)) %>% pull(Orthogroup)

intersect(imp_og, rev_og) %>% length()
length(imp_og)
length(rev_og)












